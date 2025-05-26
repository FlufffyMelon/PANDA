import click
from panda.builder import build_system


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "--config",
    "-c",
    required=True,
    type=click.Path(exists=True),
    help="Path to the config file",
)
def build(config):
    """Build an MD system from a config file."""
    build_system(config)


if __name__ == "__main__":
    cli()
