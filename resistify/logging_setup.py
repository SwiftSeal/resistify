import logging
from rich.logging import RichHandler
from rich.console import Console

console = Console()
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=console)],
)
log = logging.getLogger("rich")