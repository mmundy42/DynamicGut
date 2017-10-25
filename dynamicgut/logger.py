import logging

logger = logging.getLogger('DynamicGut')
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
handler.setFormatter(formatter)
handler.setLevel(logging.WARNING)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def log_to_file(file_name, level=logging.INFO, console=False):
    """ Turn on logging to a file.

    Parameters
    ----------
    file_name : str
        Path to file for storing log output
    level : int, optional
        Log level value (0-50)
    console : boolean, optional
        When True, keep logging to the console
    """

    fh = logging.FileHandler(file_name)
    fh.setFormatter(formatter)
    fh.setLevel(level)
    logger.addHandler(fh)
    if not console:
        ch = logger.handlers[0]
        logger.removeHandler(ch)
