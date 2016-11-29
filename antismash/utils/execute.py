"""Executable helpers for antiSMASH"""
import logging
import subprocess
try:
    from cStrinIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO
from helperlibs.wrappers.io import TemporaryFile
from antismash.config import get_config


# Ignore the pylint warning about input being redifined, as we're just
# following the subprocess names here.
# pylint: disable=redefined-builtin
def run(commands, input=None):
    """Execute commands in a system-independent manner"""

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as os_err:
        logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, os_err)
        raise
# pylint: enable=redefined-builtin


def run_hmmsearch(query_hmmfile, target_sequence, use_tempfile=False):
    """Run hmmsearch"""
    config = get_config()
    command = ["hmmsearch", "--cpu", str(config.cpus),
               query_hmmfile]
    try:
        if use_tempfile:
            with TemporaryFile(suffix='.fa') as tmp:
                command.append(tmp.name)
                with open(tmp.name, 'w') as handle:
                    handle.write(target_sequence)
                out, err, retcode = run(command)
        else:
            command.append('-')
            out, err, retcode = run(command, input=target_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmsearch returned %d: %r while searching %r', retcode,
                      err, query_hmmfile)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results


def run_hmmscan(target_hmmfile, query_sequence, opts=None):
    """Run hmmscan"""
    config = get_config()
    command = ["hmmscan", "--cpu", str(config.cpus), "--nobias"]
    if opts is not None:
        command.extend(opts)
    command.extend([target_hmmfile, '-'])
    try:
        out, err, retcode = run(command, input=query_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmscan returned %d: %r while scanning %r', retcode,
                      err, query_sequence)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results


def run_hmmpfam2(query_hmmfile, target_sequence):
    """Run hmmpfam2"""
    config = get_config()
    command = ["hmmpfam2", "--cpu", str(config.cpus),
               query_hmmfile, '-']
    try:
        out, err, retcode = run(command, input=target_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmpfam2 returned %d: %r while searching %r', retcode,
                      err, query_hmmfile)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer2-text'))
    return results


def run_hmmpress(hmmfile):
    """Run hmmpress"""
    command = ['hmmpress', hmmfile]
    try:
        out, err, retcode = run(command)
    except OSError as os_err:
        retcode = 1
        err = str(os_err)
    return out, err, retcode
