import logging

logger = logging.getLogger(__name__)


class ProgressLogger:
    def __init__(self, total_count):
        self.total_count = total_count
        self.current_count = 0
        self.last_reported_percent = -1  # Initialize with an invalid percentage

    def update(self):
        self.current_count += 1
        if self.total_count < 10:
            # For small totals, report as "n of total"
            logger.info(f"{self.current_count} of {self.total_count} complete")
        else:
            # Calculate percentage
            percent_complete = int((self.current_count / self.total_count) * 100)
            if self.current_count == self.total_count:
                logger.info("100% complete")
                self.last_reported_percent = 100
            elif (
                percent_complete % 10 == 0
                and percent_complete > 0
                and percent_complete > self.last_reported_percent
            ):
                logger.info(f"{percent_complete}% complete")
                self.last_reported_percent = percent_complete
