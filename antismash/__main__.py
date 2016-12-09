#!/usr/bin/env python
'''CLI interface for antiSMASH'''
from __future__ import print_function

import sys
import logging
import argparse
import multiprocessing

from antismash import __version__, utils
from antismash.log import setup_logging
from antismash.config import load_config, set_config
from antismash.plugins import load_plugins, check_prereqs, get_versions
from antismash.core import run_antismash


def main():
    '''Configure and parse the CLI'''
    parser = argparse.ArgumentParser()

    # we actually want at least 1 sequence file, but don't want to check here
    # this allows us to run --check-prereqs and stuff.
    parser.add_argument('sequences',
                        metavar='sequence',
                        nargs='*',
                        help="GenBank/EMBL/FASTA formated file(s)")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        default=False,
                        help="Print verbose status information to stderr.")
    parser.add_argument('-d', '--debug',
                        dest='debug',
                        action='store_true',
                        default=False,
                        help="Print debugging information to stderr.")
    parser.add_argument('--logfile',
                        dest='logfile',
                        default=argparse.SUPPRESS,
                        help="Also write logging output to a file.")
    parser.add_argument('-o', '--outputfolder',
                        dest='outputfoldername',
                        default=argparse.SUPPRESS,
                        help='Directory to write results to.')
    parser.add_argument('--cpus',
                        dest='cpus',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help='How many CPUs to use in parallel (default: {})'.format(
                            multiprocessing.cpu_count()))
    parser.add_argument('--list-plugins',
                        dest='list_plugins',
                        action='store_true',
                        default=False,
                        help="List all available sec. met. detection modules.")
    parser.add_argument('--check-prereqs',
                        dest='check_prereqs_only',
                        action='store_true',
                        default=False,
                        help="Just check if all prerequisites are met.")
    parser.add_argument('-V', '--version',
                        dest='version',
                        action='store_true',
                        default=False,
                        help="Display the version number and exit.")


    options = parser.parse_args()

    setup_logging(options)

    #if -V, show version texts and exit
    if options.version and not (options.debug or options.verbose):
        print("fake antiSMASH {}".format(__version__))
        sys.exit(0)

    load_config(options)
    set_config(options)

    #Load plugins
    logging.debug("Loading available plugins")
    plugins = load_plugins()

    if options.list_plugins:
        list_plugins(plugins)
        sys.exit(0)

    if options.version and (options.debug or options.verbose):
        versions = get_versions(plugins, options)
        print("fake antiSMASH {}".format(__version__))
        print("Identified dependencies:")
        for version in versions:
            print("    {}".format(version))
        sys.exit(0)

    #Check prerequisites
    failures = check_prereqs(plugins, options)
    if len(failures) > 0:
        for msg in failures:
            logging.error(msg)
        logging.error("Not all prerequisites met")
        sys.exit(1)
    if options.check_prereqs_only:
        logging.info("All prerequisites are met")
        sys.exit(0)

    if not options.sequences:
        parser.error("Please specify at least one sequence file")

    utils.create_outputfolder(options)

    seq_records = utils.parse_input_sequences(options)
    run_antismash(seq_records, plugins, options)


def list_plugins(plugins):
    '''List all available plugins'''
    print('Available input format plugins:')
    for plugin in plugins['input']:
        print('  * {}'.format(plugin.short_description))

    print('Available secondary metabolite cluster detection plugins:')
    for plugin in plugins['cluster']:
        print('  * {}'.format(plugin.short_description))

    print('Available cluster-specific prediction plugins:')
    for plugin in plugins['specific']:
        print('  * {}'.format(plugin.short_description))

    print('Available generic secondary metabolite prediction plugins:')
    for plugin in plugins['generic']:
        print('  * {}'.format(plugin.short_description))

    print('Available genome wide prediction plugins:')
    for plugin in plugins['genome_wide']:
        print('  * {}'.format(plugin.short_description))

    print('Available output plugins:')
    for plugin in plugins['output']:
        print('  * {}'.format(plugin.short_description))


if __name__ == '__main__':
    main()
