import exac
import typer

manager = typer.Typer()


@manager.command("hello")
def hello():
    print("hello")


@manager.command("load_db")
def load_db():
    exac.load_db()


@manager.command("load_base_coverage")
def load_base_coverage():
    exac.load_base_coverage()


@manager.command("load_variants_file")
def load_variants_file():
    exac.load_variants_file()


@manager.command("reload_variants")
def reload_variants():
    exac.load_variants_file()
    exac.load_mnps()
    exac.precalculate_metrics()


@manager.command("load_gene_models")
def load_gene_models():
    exac.load_gene_models()


@manager.command("load_cnv_models")
def load_cnv_models():
    exac.load_cnv_models()


@manager.command("load_cnv_genes")
def load_cnv_genes():
    exac.load_cnv_genes()


@manager.command("drop_cnv_genes")
def drop_cnv_genes():
    exac.drop_cnv_genes()


@manager.command("load_dbsnp_file")
def load_dbsnp_file():
    exac.load_dbsnp_file()


@manager.command("load_constraint_information")
def load_constraint_information():
    exac.load_constraint_information()


@manager.command("load_mnps")
def load_mnps():
    exac.load_mnps()


@manager.command("create_cache")
def create_cache():
    exac.create_cache()


@manager.command("precalculate_metrics")
def precalculate_metrics():
    exac.precalculate_metrics()


if __name__ == "__main__":
    manager()
