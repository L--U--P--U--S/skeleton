'''IO helpers for antiSMASH'''
import logging
import os
from os import path
from helperlibs.bio import seqio
from antismash.utils.errors import (
    NoInputFileError,
    NCBIDownloadError,
)


def parse_input_sequences(options):
    '''Parse the input sequence(s) from the given filename(s)'''
    logging.info('Parsing input sequences %r', options.sequences)

    sequences = []
    for filename in options.sequences:
        if not path.exists(filename):
            logging.error('No sequence file found at %r', filename)
            raise NoInputFileError(filename)

        with open(filename, 'r') as handle:
            content = handle.read()
            if 'Resource temporarily unavailable' in content or \
               '<h1>ServerError</h1>' in content or \
               'NCBI - WWW Error' in content:
                logging.error('NCBI server temorarily not available: downloading %s failed.',
                              path.basename(filename))
                raise NCBIDownloadError(filename)

            handle.seek(0)
            try:
                record_list = list(seqio.parse(handle))
                # TODO: filter out short sequence records
                sequences.extend(record_list)
            except ValueError as err:
                logging.error('Parsing %r failed: %s', filename, err)
                raise
            except AssertionError as err:
                logging.error('Parsing %r failed: %s', filename, err)
                raise
            except Exception as err:
                logging.error('Parsing %r failed with unhandled exception: %s', filename, err)
                raise

    # TODO: Add cleanup logic here
    return sequences


def create_outputfolder(options):
    '''create the output folder'''
    if 'outputfoldername' not in options:
        options.outputfoldername = path.splitext(path.basename(options.sequences[0]))[0]
    if not path.exists(options.outputfoldername):
        os.mkdir(options.outputfoldername)
    options.full_outputfolder_path = path.abspath(options.outputfoldername)
