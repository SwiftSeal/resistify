import subprocess
from resistify._loguru import logger


def check_dependencies() -> bool:
    try:
        from resistify.tmbed import tmbed  # noqa: F401
        from resistify.coconat import coconat  # noqa: F401
        from resistify.nlrexpress import nlrexpress  # noqa: F401
        from resistify.hmmsearch import hmmsearch  # noqa: F401
        from resistify.draw import draw  # noqa: F401
    except ImportError as e:
        logger.error(
            f"Missing dependency: {e.name}. Please install it with 'pip install resistify[full]'."
        )
        return False

    # Check if the required executables are available
    required_executables = ["jackhmmer", "hmmsearch"]
    for exe in required_executables:
        try:
            logger.info(f"Checking for {exe}...")
            subprocess.run(
                [exe, "-h"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
        except subprocess.CalledProcessError:
            logger.error(
                f"Executable '{exe}' not found. Please ensure it is installed and available in your PATH."
            )
            return False

    logger.info("All dependencies are satisfied.")
    return True
