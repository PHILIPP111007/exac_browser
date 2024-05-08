import os
import sys
import time
import gzip
import random
import itertools
import traceback
from collections import defaultdict
from multiprocessing import Process

import numpy
import pymongo
import pysam
from fastapi import Request, Response

from modules.logger import logger
from modules.settings import settings
from modules.db import get_db
from modules.templates import templates
from modules import lookups
from modules.parsing import (
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
from modules.utils import AF_BUCKETS, add_transcript_coordinate_to_variants


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
        sys.exit(0)

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
            page_content = get_gene_page_content(gene_id)
        except Exception as e:
            logger.info(e)
            continue
        f = open(os.path.join(settings.GENE_CACHE_DIR, "{}.html".format(gene_id)), "w")
        f.write(page_content)
        f.close()


def get_gene_page_content(request: Request, gene_id: str):
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
            transcript["xstart"] - settings.EXON_PADDING,
            transcript["xstop"] + settings.EXON_PADDING,
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


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if not g.db_conn:
#         g.db_conn.close()
