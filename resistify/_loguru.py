import os
import sys
from loguru import logger as _logger


class LoggerDelegator:
    def __init__(self, logger):
        self._logger = logger
        self._configure_logger()

    def _configure_logger(self):
        # Remove default handler
        self._logger.remove()

        # Add a new handler with a default format and dynamic level
        log_level = os.getenv("LOG_LEVEL", "INFO")
        self._logger.add(
            sink=sys.stdout,
            format="[{time:HH:mm:ss}] <level>{level: <8}</level> {message}",
            level=log_level,
        )

    def update_level(self, level):
        os.environ["LOG_LEVEL"] = level  # Set globally for spawned processes
        self._configure_logger()  # Reconfigure logger in the main process

    def __getattr__(self, attr):
        return getattr(self._logger, attr)


logger = LoggerDelegator(_logger)
