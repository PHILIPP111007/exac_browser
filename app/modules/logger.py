import sys
import logging

import uvicorn


def _get_logger():
    logger = logging.getLogger("uvicorn.access")
    logger.setLevel(logging.DEBUG)
    console_formatter = uvicorn.logging.ColourizedFormatter(
        "{levelprefix} {message}", style="{", use_colors=True
    )
    stream_handler = logging.StreamHandler(sys.stdout)
    logger.addHandler(stream_handler)
    logger.handlers[0].setFormatter(console_formatter)
    return logger


logger = _get_logger()
