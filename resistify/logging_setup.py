import logging
from rich.logging import RichHandler
from rich.console import Console

console = Console()
logging.basicConfig(
    level="NOTSET",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console)],
)
log = logging.getLogger("rich")
log.info("Welcome to Resistify version 0.2.0!")