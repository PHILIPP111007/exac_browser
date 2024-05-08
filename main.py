import modules.lookups as lookups
import sqlite3
import traceback

# NEW #################################################################
#######################################################################
import os
from pathlib import Path
import re
from collections import OrderedDict
import json

import uvicorn
from fastapi import FastAPI, Request, Response
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles

from modules.utils import (
    add_transcript_coordinate_to_variants,
    add_consequence_to_variant,
    get_proper_hgvs,
    remove_extraneous_vep_annotations,
    order_vep_by_csq,
    get_xpos,
)

from modules.fastapi_globals import g
from modules.logger import logger
from modules.settings import settings
from modules.db import get_db
from modules.commands import get_gene_page_content
from modules.templates import templates

#######################################################################
#######################################################################


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

static_dir = Path(__file__).parent
static_dir = static_dir.joinpath("static")
app.mount("/static", StaticFiles(directory=static_dir), name="static")


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
        except Exception as e:
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
        return get_gene_page_content(request, gene_id)


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
            transcript["xstart"] - settings.EXON_PADDING,
            transcript["xstop"] + settings.EXON_PADDING,
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
        if start is None or stop - start > settings.REGION_LIMIT or stop < start:
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
    uvicorn.run(app=settings.APP, host=settings.HOST, port=settings.PORT, reload=True)
