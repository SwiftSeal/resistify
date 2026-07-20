import os
import multiprocessing


def get_threads() -> int:
    if hasattr(os, "sched_getaffinity"):
        try:
            return len(os.sched_getaffinity(0))
        except Exception:
            pass
    return int(multiprocessing.cpu_count())
