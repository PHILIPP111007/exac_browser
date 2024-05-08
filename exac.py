import pymongo
import pysam
import lookups


from multiprocessing import Process
import sqlite3
import traceback


# NEW #################################################################
#######################################################################
import os
import sys
import glob
import time
import random
import re
import itertools
from collections import defaultdict, OrderedDict
import json
import gzip
import logging

import uvicorn
from fastapi import FastAPI, Request, Response
from fastapi.responses import HTMLResponse, RedirectResponse
from pydantic_settings import BaseSettings
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

import numpy

from parsing import (
    get_base_coverage_from_file,
    get_variants_from_sites_vcf,
    get_mnp_data,
    get_constraint_information,
    get_canonical_transcripts,
    get_omim_associations,
    get_genes_from_gencode_gtf,
    get_transcripts_from_gencode_gtf,
    get_exons_from_gencode_gtf,
    get_cnvs_from_txt,
    get_cnvs_per_gene,
    get_dbnsfp_info,
    get_snp_from_dbsnp_file,
)
from utils import (
    AF_BUCKETS,
    add_transcript_coordinate_to_variants,
    add_consequence_to_variant,
    get_proper_hgvs,
    remove_extraneous_vep_annotations,
    order_vep_by_csq,
    get_xpos,
)

from fastapi_globals import g

#######################################################################
#######################################################################

logger = logging.getLogger("uvicorn.access")
logger.setLevel(logging.DEBUG)
console_formatter = uvicorn.logging.ColourizedFormatter(
    "{asctime} {levelprefix} {message}", style="{", use_colors=True
)
stream_handler = logging.StreamHandler(sys.stdout)
logger.addHandler(stream_handler)
logger.handlers[0].setFormatter(console_formatter)


EXAC_FILES_DIRECTORY = "/Users/phil/Downloads/exac_data"
REGION_LIMIT = 1e5
EXON_PADDING = 50


class Settings(BaseSettings):
    DB_HOST: str = "localhost"
    DB_PORT: int = 27017
    DB_NAME: str = "exac"
    DEBUG: bool = True
    SECRET_KEY: str = "development key"
    LOAD_DB_PARALLEL_PROCESSES: int = 4  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    SITES_VCFS: list[str] = glob.glob(
        os.path.join(EXAC_FILES_DIRECTORY, "ExAC*.vcf.gz")
    )
    GENCODE_GTF: str = os.path.join(EXAC_FILES_DIRECTORY, "gencode.gtf.gz")
    CANONICAL_TRANSCRIPT_FILE: str = os.path.join(
        EXAC_FILES_DIRECTORY,
        "canonical_transcripts.txt.gz",
    )
    OMIM_FILE: str = os.path.join(EXAC_FILES_DIRECTORY, "omim_info.txt.gz")
    BASE_COVERAGE_FILES: list[str] = glob.glob(
        os.path.join(
            EXAC_FILES_DIRECTORY,
            "coverage",
            "Panel.*.coverage.txt.gz",
        )
    )
    DBNSFP_FILE: str = os.path.join(EXAC_FILES_DIRECTORY, "dbNSFP2.6_gene.gz")
    CONSTRAINT_FILE: str = os.path.join(
        EXAC_FILES_DIRECTORY,
        "forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz",
    )
    MNP_FILE: str = os.path.join(
        EXAC_FILES_DIRECTORY,
        "MNPs_NotFiltered_ForBrowserRelease.txt.gz",
    )
    CNV_FILE: str = os.path.join(
        EXAC_FILES_DIRECTORY,
        "exac-gencode-exon.cnt.final.pop3",
    )
    CNV_GENE_FILE: str = os.path.join(EXAC_FILES_DIRECTORY, "exac-final-cnvs.gene.rank")
    # How to get a dbsnp142.txt.bgz file:
    #   wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/database/organism_data/b142_SNPChrPosOnRef_105.bcp.gz
    #   zcat b142_SNPChrPosOnRef_105.bcp.gz | awk '$3 != ""' | perl -pi -e 's/ +/\t/g' | sort -k2,2 -k3,3n | bgzip -c > dbsnp142.txt.bgz
    #   tabix -s 2 -b 3 -e 3 dbsnp142.txt.bgz
    DBSNP_FILE: str = os.path.join(EXAC_FILES_DIRECTORY, "dbsnp142.txt.bgz")
    READ_VIZ_DIR: str = os.path.join("/mongo", "readviz")

    GENE_CACHE_DIR: str = os.path.join(os.path.dirname(__file__), "gene_cache")
    GENES_TO_CACHE: set[str] = {
        l.strip("\n")
        for l in open(os.path.join(os.path.dirname(__file__), "genes_to_cache.txt"))
    }


settings = Settings()


app = FastAPI(
    title="EXAC Browser",
    description="Description",
    version="0.0.1",
    terms_of_service=None,
    contact=None,
    license_info=None,
    docs_url="/api/docs/",
    redoc_url="/api/redoc/",
    openapi_url="/api/openapi.json",
)
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")


def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=settings.DB_HOST, port=settings.DB_PORT)
    return client[settings.DB_NAME]


def parse_tabix_file_subset(tabix_filenames, subset_i, subset_n, record_parser):
    """
    Returns a generator of parsed record objects (as returned by record_parser) for the i'th out n subset of records
    across all the given tabix_file(s). The records are split by files and contigs within files, with 1/n of all contigs
    from all files being assigned to this the i'th subset.

    Args:
        tabix_filenames: a list of one or more tabix-indexed files. These will be opened using pysam.Tabixfile
        subset_i: zero-based number
        subset_n: total number of subsets
        record_parser: a function that takes a file-like object and returns a generator of parsed records
    """
    start_time = time.time()
    open_tabix_files = [
        pysam.Tabixfile(tabix_filename) for tabix_filename in tabix_filenames
    ]
    tabix_file_contig_pairs = [
        (tabix_file, contig)
        for tabix_file in open_tabix_files
        for contig in tabix_file.contigs
    ]
    tabix_file_contig_subset = tabix_file_contig_pairs[
        subset_i::subset_n
    ]  # get every n'th tabix_file/contig pair
    short_filenames = ", ".join(map(os.path.basename, tabix_filenames))
    num_file_contig_pairs = len(tabix_file_contig_subset)
    logger.info(
        (
            "Loading subset %(subset_i)s of %(subset_n)s total: %(num_file_contig_pairs)s contigs from "
            "%(short_filenames)s"
        )
        % locals()
    )
    counter = 0
    for tabix_file, contig in tabix_file_contig_subset:
        header_iterator = tabix_file.header
        records_iterator = tabix_file.fetch(contig, 0, 10**9, multiple_iterators=True)
        for parsed_record in record_parser(
            itertools.chain(header_iterator, records_iterator)
        ):
            counter += 1
            yield parsed_record

            if counter % 100000 == 0:
                seconds_elapsed = int(time.time() - start_time)
                logger.info(
                    (
                        "Loaded %(counter)s records from subset %(subset_i)s of %(subset_n)s from %(short_filenames)s "
                        "(%(seconds_elapsed)s seconds)"
                    )
                    % locals()
                )

    logger.info(
        "Finished loading subset %(subset_i)s from  %(short_filenames)s (%(counter)s records)"
        % locals()
    )


def load_base_coverage():
    def load_coverage(coverage_files, i, n, db):
        coverage_generator = parse_tabix_file_subset(
            coverage_files, i, n, get_base_coverage_from_file
        )
        try:
            db.base_coverage.insert(coverage_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when coverage_generator is empty

    db = get_db()
    db.base_coverage.drop()
    logger.info("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    db.base_coverage.create_index("xpos")

    procs = []
    coverage_files = settings.BASE_COVERAGE_FILES
    num_procs = settings.LOAD_DB_PARALLEL_PROCESSES
    random.shuffle(settings.BASE_COVERAGE_FILES)
    for i in range(num_procs):
        p = Process(target=load_coverage, args=(coverage_files, i, num_procs, db))
        p.start()
        procs.append(p)
    return procs

    # logger.info 'Done loading coverage. Took %s seconds' % int(time.time() - start_time)


def load_variants_file():
    def load_variants(sites_file, i, n, db):
        variants_generator = parse_tabix_file_subset(
            [sites_file], i, n, get_variants_from_sites_vcf
        )
        try:
            db.variants.insert(variants_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when variant_generator is empty

    db = get_db()
    db.variants.drop()
    logger.info("Dropped db.variants")

    # grab variants from sites VCF
    db.variants.create_index("xpos")
    db.variants.create_index("xstart")
    db.variants.create_index("xstop")
    db.variants.create_index("rsid")
    db.variants.create_index("genes")
    db.variants.create_index("transcripts")

    sites_vcfs = settings.SITES_VCFS
    if len(sites_vcfs) == 0:
        raise IOError("No vcf file found")
    elif len(sites_vcfs) > 1:
        raise Exception("More than one sites vcf file found: %s" % sites_vcfs)

    num_procs = settings.LOAD_DB_PARALLEL_PROCESSES
    load_variants(sites_vcfs[0], 0, num_procs, db)


def load_constraint_information():
    db = get_db()

    db.constraint.drop()
    logger.info("Dropped db.constraint.")

    start_time = time.time()

    with gzip.open(settings.CONSTRAINT_FILE) as constraint_file:
        for transcript in get_constraint_information(constraint_file):
            db.constraint.insert(transcript, w=0)

    db.constraint.create_index("transcript")
    logger.info(
        "Done loading constraint info. Took %s seconds" % int(time.time() - start_time)
    )


def load_mnps():
    db = get_db()
    start_time = time.time()

    db.variants.create_index("has_mnp")
    logger.info("Done indexing.")
    while db.variants.find_and_modify(
        {"has_mnp": True}, {"$unset": {"has_mnp": "", "mnps": ""}}
    ):
        pass
    logger.info("Deleted MNP data.")

    with gzip.open(settings.MNP_FILE) as mnp_file:
        for mnp in get_mnp_data(mnp_file):
            variant = lookups.get_raw_variant(
                db, mnp["xpos"], mnp["ref"], mnp["alt"], True
            )
            db.variants.find_and_modify(
                {"_id": variant["_id"]},
                {"$set": {"has_mnp": True}, "$push": {"mnps": mnp}},
                w=0,
            )

    db.variants.create_index("has_mnp")
    logger.info(
        "Done loading MNP info. Took %s seconds" % int(time.time() - start_time)
    )


def load_gene_models():
    db = get_db()

    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    logger.info("Dropped db.genes, db.transcripts, and db.exons.")

    start_time = time.time()

    canonical_transcripts = {}
    with gzip.open(settings.CANONICAL_TRANSCRIPT_FILE) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript

    omim_annotations = {}
    with gzip.open(settings.OMIM_FILE) as omim_file:
        for fields in get_omim_associations(omim_file):
            logger.info(fields)
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)

    dbnsfp_info = {}
    with gzip.open(settings.DBNSFP_FILE) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [
                other_name.upper() for other_name in dbnsfp_gene["gene_other_names"]
            ]
            dbnsfp_info[dbnsfp_gene["ensembl_gene"]] = (
                dbnsfp_gene["gene_full_name"],
                other_names,
            )

    logger.info(
        "Done loading metadata. Took %s seconds" % int(time.time() - start_time)
    )

    # grab genes from GTF
    start_time = time.time()
    with gzip.open(settings.GENCODE_GTF) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene["gene_id"]
            if gene_id in canonical_transcripts:
                gene["canonical_transcript"] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene["omim_accession"] = omim_annotations[gene_id][0]
                gene["omim_description"] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene["full_gene_name"] = dbnsfp_info[gene_id][0]
                gene["other_names"] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)

    logger.info("Done loading genes. Took %s seconds" % int(time.time() - start_time))

    start_time = time.time()
    db.genes.create_index("gene_id")
    db.genes.create_index("gene_name_upper")
    db.genes.create_index("gene_name")
    db.genes.create_index("other_names")
    db.genes.create_index("xstart")
    db.genes.create_index("xstop")
    logger.info(
        "Done indexing gene table. Took %s seconds" % int(time.time() - start_time)
    )

    # and now transcripts
    start_time = time.time()
    with gzip.open(settings.GENCODE_GTF) as gtf_file:
        db.transcripts.insert(
            (transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)),
            w=0,
        )
    logger.info(
        "Done loading transcripts. Took %s seconds" % int(time.time() - start_time)
    )

    start_time = time.time()
    db.transcripts.create_index("transcript_id")
    db.transcripts.create_index("gene_id")
    logger.info(
        "Done indexing transcript table. Took %s seconds"
        % int(time.time() - start_time)
    )

    # Building up gene definitions
    start_time = time.time()
    with gzip.open(settings.GENCODE_GTF) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    logger.info("Done loading exons. Took %s seconds" % int(time.time() - start_time))

    start_time = time.time()
    db.exons.create_index("exon_id")
    db.exons.create_index("transcript_id")
    db.exons.create_index("gene_id")
    logger.info(
        "Done indexing exon table. Took %s seconds" % int(time.time() - start_time)
    )

    return []


def load_cnv_models():
    db = get_db()

    db.cnvs.drop()
    logger.info("Dropped db.cnvs.")

    start_time = time.time()
    with open(settings.CNV_FILE) as cnv_txt_file:
        for cnv in get_cnvs_from_txt(cnv_txt_file):
            db.cnvs.insert(cnv, w=0)
            # progress.update(gtf_file.fileobj.tell())
        # progress.finish()

    logger.info("Done loading CNVs. Took %s seconds" % int(time.time() - start_time))


def drop_cnv_genes():
    db = get_db()
    start_time = time.time()
    db.cnvgenes.drop()


def load_cnv_genes():
    db = get_db()
    start_time = time.time()
    with open(settings.CNV_GENE_FILE) as cnv_gene_file:
        for cnvgene in get_cnvs_per_gene(cnv_gene_file):
            db.cnvgenes.insert(cnvgene, w=0)
            # progress.update(gtf_file.fileobj.tell())
        # progress.finish()

    logger.info(
        "Done loading CNVs in genes. Took %s seconds" % int(time.time() - start_time)
    )


def load_dbsnp_file():
    db = get_db()

    def load_dbsnp(dbsnp_file, i, n, db):
        if os.path.isfile(dbsnp_file + ".tbi"):
            dbsnp_record_generator = parse_tabix_file_subset(
                [dbsnp_file], i, n, get_snp_from_dbsnp_file
            )
            try:
                db.dbsnp.insert(dbsnp_record_generator, w=0)
            except pymongo.errors.InvalidOperation:
                pass  # handle error when coverage_generator is empty

        else:
            with gzip.open(dbsnp_file) as f:
                db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(f)), w=0)

    db.dbsnp.drop()
    db.dbsnp.create_index("rsid")
    db.dbsnp.create_index("xpos")
    start_time = time.time()
    dbsnp_file = settings.DBSNP_FILE

    logger.info("Loading dbsnp from %s" % dbsnp_file)
    if os.path.isfile(dbsnp_file + ".tbi"):
        num_procs = settings.LOAD_DB_PARALLEL_PROCESSES
    else:
        # see if non-tabixed .gz version exists
        if os.path.isfile(dbsnp_file):
            logger.info(
                (
                    "WARNING: %(dbsnp_file)s.tbi index file not found. Will use single thread to load dbsnp."
                    "To create a tabix-indexed dbsnp file based on UCSC dbsnp, do: \n"
                    "   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz \n"
                    "   gzcat snp141.txt.gz | cut -f 1-5 | bgzip -c > snp141.txt.bgz \n"
                    "   tabix -0 -s 2 -b 3 -e 4 snp141.txt.bgz"
                )
                % locals()
            )
            num_procs = 1
        else:
            raise Exception("dbsnp file %s(dbsnp_file)s not found." % locals())

    procs = []
    for i in range(num_procs):
        p = Process(target=load_dbsnp, args=(dbsnp_file, i, num_procs, db))
        p.start()
        procs.append(p)

    return procs
    # logger.info 'Done loading dbSNP. Took %s seconds' % int(time.time() - start_time)

    # start_time = time.time()
    # db.dbsnp.create_index('rsid')
    # logger.info 'Done indexing dbSNP table. Took %s seconds' % int(time.time() - start_time)


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = input(
        "This will drop the database and reload. Are you sure you want to continue? [y/n] "
    )
    if confirm != "y":
        logger.info("Exiting...")
        sys.exit(1)

    for load_function in [
        load_variants_file,
        load_dbsnp_file,
        load_base_coverage,
        load_gene_models,
        load_constraint_information,
        load_cnv_models,
        load_cnv_genes,
    ]:
        load_function()

    logger.info("Done! Loading MNPs...")
    load_mnps()
    logger.info("Done! Creating cache...")
    create_cache()
    logger.info("Done!")


def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    for gene in get_db().genes.find():
        autocomplete_strings.append(gene["gene_name"])
        if "other_names" in gene:
            autocomplete_strings.extend(gene["other_names"])
    f = open(os.path.join(os.path.dirname(__file__), "autocomplete_strings.txt"), "w")
    for s in sorted(autocomplete_strings):
        f.write(s + "\n")
    f.close()

    # create static gene pages for genes in
    if not os.path.exists(settings.GENE_CACHE_DIR):
        os.makedirs(settings.GENE_CACHE_DIR)

    # get list of genes ordered by num_variants
    for gene_id in settings.GENES_TO_CACHE:
        try:
            page_content = _get_gene_page_content(gene_id)
        except Exception as e:
            logger.info(e)
            continue
        f = open(os.path.join(settings.GENE_CACHE_DIR, "{}.html".format(gene_id)), "w")
        f.write(page_content)
        f.close()


def precalculate_metrics():
    db = get_db()
    logger.info("Reading %s variants..." % db.variants.estimated_document_count())
    metrics = defaultdict(list)
    binned_metrics = defaultdict(list)
    progress = 0
    start_time = time.time()
    for variant in db.variants.find():
        for metric, value in variant["quality_metrics"].iteritems():
            metrics[metric].append(float(value))
        qual = float(variant["site_quality"])
        metrics["site_quality"].append(qual)
        if variant["allele_num"] == 0:
            continue
        if variant["allele_count"] == 1:
            binned_metrics["singleton"].append(qual)
        elif variant["allele_count"] == 2:
            binned_metrics["doubleton"].append(qual)
        else:
            for af in AF_BUCKETS:
                if float(variant["allele_count"]) / variant["allele_num"] < af:
                    binned_metrics[af].append(qual)
                    break
        progress += 1
        if not progress % 100000:
            logger.info(
                "Read %s variants. Took %s seconds"
                % (progress, int(time.time() - start_time))
            )
    logger.info("Done reading variants. Dropping metrics database... ")
    db.metrics.drop()
    logger.info("Dropped metrics database. Calculating metrics...")
    for metric in metrics:
        bin_range = None
        data = map(numpy.log, metrics[metric]) if metric == "DP" else metrics[metric]
        if metric == "FS":
            bin_range = (0, 20)
        elif metric == "VQSLOD":
            bin_range = (-20, 20)
        elif metric == "InbreedingCoeff":
            bin_range = (0, 1)
        if bin_range is not None:
            data = [x if (x > bin_range[0]) else bin_range[0] for x in data]
            data = [x if (x < bin_range[1]) else bin_range[1] for x in data]
        hist = numpy.histogram(data, bins=40, range=bin_range)
        edges = hist[1]
        # mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
        lefts = [edges[i] for i in range(len(edges) - 1)]
        db.metrics.insert({"metric": metric, "mids": lefts, "hist": list(hist[0])})
    for metric in binned_metrics:
        hist = numpy.histogram(map(numpy.log, binned_metrics[metric]), bins=40)
        edges = hist[1]
        mids = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
        db.metrics.insert(
            {"metric": "binned_%s" % metric, "mids": mids, "hist": list(hist[0])}
        )
    db.metrics.create_index("metric")
    logger.info("Done pre-calculating metrics!")


def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    g.db_conn = connect_db()
    return g.db_conn


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if not g.db_conn:
#         g.db_conn.close()


@app.get("/", response_class=HTMLResponse)
def homepage(request: Request):
    return templates.TemplateResponse(request=request, name="homepage.html")


@app.get("/autocomplete/<query>")
def awesome_autocomplete(request: Request, query):
    if not hasattr(g, "autocomplete_strings"):
        g.autocomplete_strings = [
            s.strip()
            for s in open(
                os.path.join(os.path.dirname(__file__), "autocomplete_strings.txt")
            )
        ]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(
        json.dumps([{"value": s} for s in suggestions]), media_type="application/json"
    )


@app.get("/awesome/")
def awesome(request: Request, query: str):
    db = get_db()
    datatype, identifier = lookups.get_awesomebar_result(db, query)

    logger.info("Searched for %s: %s" % (datatype, identifier))
    if datatype == "gene":
        return RedirectResponse("/gene/{}".format(identifier))
    elif datatype == "transcript":
        return RedirectResponse("/transcript/{}".format(identifier))
    elif datatype == "variant":
        return RedirectResponse("/variant/{}".format(identifier))
    elif datatype == "region":
        return RedirectResponse("/region/{}".format(identifier))
    elif datatype == "dbsnp_variant_set":
        return RedirectResponse("/dbsnp/{}".format(identifier))
    elif datatype == "error":
        return RedirectResponse("/error/{}".format(identifier))
    elif datatype == "not_found":
        return RedirectResponse("/not_found/{}".format(identifier))
    else:
        raise Exception


@app.get("/variant/<variant_str>/", response_class=HTMLResponse)
def variant_page(request: Request, variant_str: str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split("-")
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups.get_variant(db, xpos, ref, alt)

        if variant is None:
            variant = {"chrom": chrom, "pos": pos, "xpos": xpos, "ref": ref, "alt": alt}
        consequences = OrderedDict()
        if "vep_annotations" in variant:
            add_consequence_to_variant(variant)
            variant["vep_annotations"] = remove_extraneous_vep_annotations(
                variant["vep_annotations"]
            )
            variant["vep_annotations"] = order_vep_by_csq(
                variant["vep_annotations"]
            )  # Adds major_consequence
            for annotation in variant["vep_annotations"]:
                annotation["HGVS"] = get_proper_hgvs(annotation)
                consequences.setdefault(annotation["major_consequence"], {}).setdefault(
                    annotation["Gene"], []
                ).append(annotation)
        base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
        any_covered = any([x["has_coverage"] for x in base_coverage])
        metrics = lookups.get_metrics(db, variant)

        # check the appropriate sqlite db to get the *expected* number of
        # available bams and *actual* number of available bams for this variant
        sqlite_db_path = os.path.join(
            settings.READ_VIZ_DIR,
            "combined_bams",
            chrom,
            "combined_chr%s_%03d.db" % (chrom, pos % 1000),
        )
        logger.info(sqlite_db_path)
        try:
            read_viz_db = sqlite3.connect(sqlite_db_path)
            n_het = read_viz_db.execute(
                "select n_expected_samples, n_available_samples from t "
                "where chrom=? and pos=? and ref=? and alt=? and het_or_hom_or_hemi=?",
                (chrom, pos, ref, alt, "het"),
            ).fetchone()
            n_hom = read_viz_db.execute(
                "select n_expected_samples, n_available_samples from t "
                "where chrom=? and pos=? and ref=? and alt=? and het_or_hom_or_hemi=?",
                (chrom, pos, ref, alt, "hom"),
            ).fetchone()
            n_hemi = read_viz_db.execute(
                "select n_expected_samples, n_available_samples from t "
                "where chrom=? and pos=? and ref=? and alt=? and het_or_hom_or_hemi=?",
                (chrom, pos, ref, alt, "hemi"),
            ).fetchone()
            read_viz_db.close()
        except Exception:
            logger.error("Error when accessing sqlite db: %s - %s", sqlite_db_path, e)
            n_het = n_hom = n_hemi = None

        read_viz_dict = {
            "het": {
                "n_expected": n_het[0]
                if n_het is not None and n_het[0] is not None
                else 0,
                "n_available": n_het[1]
                if n_het is not None and n_het[1] is not None
                else 0,
            },
            "hom": {
                "n_expected": n_hom[0]
                if n_hom is not None and n_hom[0] is not None
                else 0,
                "n_available": n_hom[1]
                if n_hom is not None and n_hom[1] is not None
                else 0,
            },
            "hemi": {
                "n_expected": n_hemi[0]
                if n_hemi is not None and n_hemi[0] is not None
                else 0,
                "n_available": n_hemi[1]
                if n_hemi is not None and n_hemi[1] is not None
                else 0,
            },
        }

        total_available = 0
        total_expected = 0
        for het_or_hom_or_hemi in ("het", "hom", "hemi"):
            total_available += read_viz_dict[het_or_hom_or_hemi]["n_available"]
            total_expected += read_viz_dict[het_or_hom_or_hemi]["n_expected"]

            read_viz_dict[het_or_hom_or_hemi]["readgroups"] = [
                "%(chrom)s-%(pos)s-%(ref)s-%(alt)s_%(het_or_hom_or_hemi)s%(i)s"
                % locals()
                for _ in range(read_viz_dict[het_or_hom_or_hemi]["n_available"])
            ]  # eg. '1-157768000-G-C_hom1',

            read_viz_dict[het_or_hom_or_hemi]["urls"] = [
                # os.path.join('combined_bams', chrom, 'combined_chr%s_%03d.bam' % (chrom, pos % 1000))
                os.path.join(
                    "combined_bams",
                    chrom,
                    "combined_chr%s_%03d.bam" % (chrom, pos % 1000),
                )
                for _ in range(read_viz_dict[het_or_hom_or_hemi]["n_available"])
            ]

        read_viz_dict["total_available"] = total_available
        read_viz_dict["total_expected"] = total_expected

        logger.info("Rendering variant: %s" % variant_str)
        return templates.TemplateResponse(
            request=request,
            name="variant.html",
            context=dict(
                variant=variant,
                base_coverage=base_coverage,
                consequences=consequences,
                any_covered=any_covered,
                metrics=metrics,
                read_viz=read_viz_dict,
            ),
        )
    except Exception:
        logger.info(
            "Failed on variant:", variant_str, ";Error=", traceback.format_exc()
        )
        return Response(status_code=404)


@app.get("/gene/<gene_id>/")
def gene_page(request: Request, gene_id: str):
    if gene_id in settings.GENES_TO_CACHE:
        return open(
            os.path.join(settings.GENE_CACHE_DIR, "{}.html".format(gene_id))
        ).read()
    else:
        return _get_gene_page_content(request, gene_id)


def _get_gene_page_content(request: Request, gene_id: str):
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if gene is None:
            return Response(status_code=404)

        variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
        transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)

        # Get some canonical transcript and corresponding info
        transcript_id = gene["canonical_transcript"]
        transcript = lookups.get_transcript(db, transcript_id)
        variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
        cnvs_in_transcript = lookups.get_exons_cnvs(db, transcript_id)
        cnvs_per_gene = lookups.get_cnvs(db, gene_id)
        coverage_stats = lookups.get_coverage_for_transcript(
            db,
            transcript["xstart"] - EXON_PADDING,
            transcript["xstop"] + EXON_PADDING,
        )
        add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
        constraint_info = lookups.get_constraint_for_transcript(db, transcript_id)

        t = templates.TemplateResponse(
            request=request,
            name="gene.html",
            context=dict(
                gene=gene,
                transcript=transcript,
                variants_in_gene=variants_in_gene,
                variants_in_transcript=variants_in_transcript,
                transcripts_in_gene=transcripts_in_gene,
                coverage_stats=coverage_stats,
                cnvs=cnvs_in_transcript,
                cnvgenes=cnvs_per_gene,
                constraint=constraint_info,
            ),
        )
        logger.info("Rendering gene: %s" % gene_id)
        return t
    except Exception:
        logger.info("Failed on gene:", gene_id, ";Error=", traceback.format_exc())
        return Response(status_code=404)


@app.get("/transcript/<transcript_id>/")
def transcript_page(request: Request, transcript_id):
    db = get_db()
    try:
        transcript = lookups.get_transcript(db, transcript_id)

        gene = lookups.get_gene(db, transcript["gene_id"])
        gene["transcripts"] = lookups.get_transcripts_in_gene(db, transcript["gene_id"])
        variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
        cnvs_in_transcript = lookups.get_exons_cnvs(db, transcript_id)
        cnvs_per_gene = lookups.get_cnvs(db, transcript["gene_id"])
        coverage_stats = lookups.get_coverage_for_transcript(
            db,
            transcript["xstart"] - EXON_PADDING,
            transcript["xstop"] + EXON_PADDING,
        )

        add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

        t = templates.TemplateResponse(
            request=request,
            name="transcript.html",
            context=dict(
                transcript=transcript,
                transcript_json=json.dumps(transcript),
                variants_in_transcript=variants_in_transcript,
                variants_in_transcript_json=json.dumps(variants_in_transcript),
                coverage_stats=coverage_stats,
                coverage_stats_json=json.dumps(coverage_stats),
                gene=gene,
                gene_json=json.dumps(gene),
                cnvs=cnvs_in_transcript,
                cnvs_json=json.dumps(cnvs_in_transcript),
                cnvgenes=cnvs_per_gene,
                cnvgenes_json=json.dumps(cnvs_per_gene),
            ),
        )

        logger.info("Rendering transcript: %s" % transcript_id)
        return t
    except Exception:
        logger.info(
            "Failed on transcript:", transcript_id, ";Error=", traceback.format_exc()
        )
        return Response(status_code=404)


@app.get("/region/<region_id>/")
def region_page(request: Request, region_id):
    db = get_db()
    try:
        region = region_id.split("-")

        chrom = region[0]
        start = None
        stop = None
        if len(region) == 3:
            chrom, start, stop = region
            start = int(start)
            stop = int(stop)
        if start is None or stop - start > REGION_LIMIT or stop < start:
            return templates.TemplateResponse(
                request=request,
                name="region.html",
                context=dict(
                    genes_in_region=None,
                    variants_in_region=None,
                    chrom=chrom,
                    start=start,
                    stop=stop,
                    coverage=None,
                ),
            )
        if start == stop:
            start -= 20
            stop += 20
        genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
        variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
        xstart = get_xpos(chrom, start)
        xstop = get_xpos(chrom, stop)
        coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
        t = templates.TemplateResponse(
            request=request,
            name="region.html",
            context=dict(
                genes_in_region=genes_in_region,
                variants_in_region=variants_in_region,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=coverage_array,
            ),
        )
        logger.info("Rendering region: %s" % region_id)
        return t
    except Exception:
        logger.info("Failed on region:", region_id, ";Error=", traceback.format_exc())
        return Response(status_code=404)


@app.get("/dbsnp/<rsid>/")
def dbsnp_page(request: Request, rsid):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        logger.info("Rendering rsid: %s" % rsid)
        return templates.TemplateResponse(
            request=request,
            name="region.html",
            context=dict(
                rsid=rsid,
                variants_in_region=variants,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=None,
                genes_in_region=None,
            ),
        )
    except Exception:
        logger.info("Failed on rsid:", rsid, ";Error=", traceback.format_exc())
        return Response(status_code=404)


@app.get("/not_found/<query>/")
def not_found_page(request: Request, query):
    return templates.TemplateResponse(
        request=request,
        name="not_found.html",
        status_code=404,
        context=dict(query=query),
    )


@app.get("/error/<query>/")
def error_page(request: Request, query):
    return templates.TemplateResponse(
        request=request, name="error.html", status_code=404, context=dict(query=query)
    )


@app.get("/downloads/")
def downloads_page(request: Request):
    return templates.TemplateResponse(request=request, name="downloads.html")


@app.get("/about/")
def about_page(request: Request):
    return templates.TemplateResponse(request=request, name="about.html")


@app.get("/participants/")
def participants_page(request: Request):
    return templates.TemplateResponse(request=request, name="about.html")


@app.get("/terms/")
def terms_page(request: Request):
    return templates.TemplateResponse(request=request, name="terms.html")


@app.get("/contact/")
def contact_page(request: Request):
    return templates.TemplateResponse(request=request, name="contact.html")


@app.get("/faq/")
def faq_page(request: Request):
    return templates.TemplateResponse(request=request, name="faq.html")


@app.get("/text/")
def text_page(request: Request, query: str):
    db = get_db()
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype in ["gene", "transcript"]:
        gene = lookups.get_gene(db, identifier)
        link = (
            "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr%(chrom)s%%3A%(start)s-%(stop)s"
            % gene
        )
        output = """Searched for %s. Found %s.
%s; Canonical: %s.
%s""" % (query, identifier, gene["full_gene_name"], gene["canonical_transcript"], link)
        output += (
            ""
            if "omim_accession" not in gene
            else """
In OMIM: %(omim_description)s
http://omim.org/entry/%(omim_accession)s"""
            % gene
        )
        return output
    elif datatype == "error" or datatype == "not_found":
        return "Gene/transcript %s not found" % query
    else:
        return "Search types other than gene transcript not yet supported"


@app.get("/read_viz/<path:path>/")
def read_viz_files(request: Request, path):
    full_path = os.path.abspath(os.path.join(settings.READ_VIZ_DIR, path))

    # security check - only files under READ_VIZ_DIR should be accsessible
    if not full_path.startswith(settings.READ_VIZ_DIR):
        return "Invalid path: %s" % path

    # handle igv.js Range header which it uses to request a subset of a .bam
    range_header = request.headers.get("Range", None)
    if not range_header:
        return send_from_directory(settings.READ_VIZ_DIR, path)

    m = re.search(r"(\d+)-(\d*)", range_header)
    if not m:
        error_msg = "ERROR: unexpected range header syntax: %s" % range_header
        logger.error(error_msg)
        return error_msg

    size = os.path.getsize(full_path)
    offset = int(m.group(1))
    length = int(m.group(2) or size) - offset

    data = None
    with open(full_path, "rb") as f:
        f.seek(offset)
        data = f.read(length)

    rv = Response(
        data, 206, media_type="application/octet-stream", direct_passthrough=True
    )
    rv.headers.add(
        "Content-Range", "bytes {0}-{1}/{2}".format(offset, offset + length - 1, size)
    )

    logger.info(
        "readviz: range request: %s-%s %s" % (m.group(1), m.group(2), full_path)
    )
    return rv


if __name__ == "__main__":
    uvicorn.run("exac:app", host="0.0.0.0", port=8000, reload=True)
