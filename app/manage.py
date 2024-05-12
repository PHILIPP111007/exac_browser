import typer

from modules.logger import logger
from modules import commands

manager = typer.Typer()


@manager.command("hello")
def hello():
    logger.info("hello")


@manager.command("load_db")
def load_db():
    commands.load_db()


@manager.command("load_base_coverage")
def load_base_coverage():
    commands.load_base_coverage()


@manager.command("load_variants_file")
def load_variants_file():
    commands.load_variants_file()


@manager.command("reload_variants")
def reload_variants():
    commands.load_variants_file()
    commands.load_mnps()
    commands.precalculate_metrics()


@manager.command("load_gene_models")
def load_gene_models():
    commands.load_gene_models()


@manager.command("load_cnv_models")
def load_cnv_models():
    commands.load_cnv_models()


@manager.command("load_cnv_genes")
def load_cnv_genes():
    commands.load_cnv_genes()


@manager.command("drop_cnv_genes")
def drop_cnv_genes():
    commands.drop_cnv_genes()


@manager.command("load_dbsnp_file")
def load_dbsnp_file():
    commands.load_dbsnp_file()


@manager.command("load_constraint_information")
def load_constraint_information():
    commands.load_constraint_information()


@manager.command("load_mnps")
def load_mnps():
    commands.load_mnps()


@manager.command("create_cache")
def create_cache():
    # commands.create_cache()
    logger.warning("now is not working because of template rendering")


@manager.command("precalculate_metrics")
def precalculate_metrics():
    commands.precalculate_metrics()


if __name__ == "__main__":
    manager()
