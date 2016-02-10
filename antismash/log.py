import os
import logging
from os import path


def setup_logging(options):
    "Set up the logging output"
    if options.debug:
        log_level = logging.DEBUG
    elif options.verbose:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    logging.basicConfig(format='%(levelname)s: %(message)s',
                        level=log_level)
    if 'logfile' in options:
        if not (path.dirname(options.logfile) == "" or os.path.exists(path.dirname(options.logfile))):
            os.mkdir(path.dirname(options.logfile))
        fh = logging.FileHandler(options.logfile)
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logging.getLogger('').addHandler(fh)
