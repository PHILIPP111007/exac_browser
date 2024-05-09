import os
import glob
from multiprocessing import cpu_count

from pydantic_settings import BaseSettings


EXAC_FILES_DIRECTORY = "/Users/phil/Downloads/exac_data"
CPU_COUNT = cpu_count()


class Settings(BaseSettings):
    APP: str = os.environ.get("APP", "main:app")
    HOST: str = os.environ.get("HOST", "0.0.0.0")
    PORT: int = os.environ.get("PORT", 8000)

    SECRET_KEY: str = os.environ.get("SECRET_KEY", "12345")
    DEBUG: bool = True

    DB_HOST: str = os.environ.get("DB_HOST", "localhost")
    DB_PORT: int = os.environ.get("DB_PORT", 27017)
    DB_NAME: str = os.environ.get("DB_NAME", "exac")

    SEARCH_LIMIT: int = 10000

    LOAD_DB_PARALLEL_PROCESSES: int = CPU_COUNT  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    SITES_VCFS: list[str] = glob.glob(
        os.path.join(EXAC_FILES_DIRECTORY, "ExAC*.vcf")
    )  # NOT GZIP
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
    REGION_LIMIT: float = 1e5
    EXON_PADDING: int = 50


settings = Settings()
