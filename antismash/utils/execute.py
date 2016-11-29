"""Executable helpers for antiSMASH"""
import logging
import subprocess
try:
    from cStrinIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings
import os
import signal
from multiprocessing import Pool

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


def run_p(options, args):
    """Execute commands in a parallel manner"""
    os.setpgid(0, 0)
    p = Pool(options.cpus)
    jobs = p.map_async(child_process, args)

    try:
        # 0xFFFF is just "a lot of seconds"
        errors = jobs.get(0xFFFF)
    except KeyboardInterrupt:
        logging.error("Interrupted by user")
        os.killpg(os.getpid(), signal.SIGTERM)
        return 128

    return errors


def child_process(args):
    """Called by multiprocessing's map or map_async method"""
    try:
        logging.debug("Calling {!r}".format(" ".join(args)))
        p = subprocess.Popen(args) # actually only a chunk of args
        p.communicate()
    except KeyboardInterrupt:
        logging.debug("{!r}: Ignoring interrupt".format(p))

    return p.returncode


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


def run_meme(meme_dir, options):
    """Set paths, check existing files and run MEME in parallel on each promoter set"""
    args = []
    for plus_minus in os.listdir(meme_dir):
        input_file = os.path.join(meme_dir, plus_minus, "promoters.fasta")
        output_file = os.path.join(meme_dir, plus_minus, "meme.xml")

        # input file present and size not zero --> should be fine to use
        if os.path.isfile(input_file) and os.path.getsize(input_file) > 0:
            # output file already present and size not zero --> do not run MEME again on this one
            if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
                pass
            else:
                args.append([
                    "meme", input_file,
                    "-oc", os.path.dirname(input_file),
                    "-dna",
                    "-nostatus",
                    "-mod", "anr",
                    "-nmotifs", "1",
                    "-minw", "6",
                    "-maxw", "12",
                    "-revcomp",
                    "-evt", "1.0e+005",
                ])

    errors = run_p(options, args)
    return sum(errors)


def run_fimo(meme_dir, fimo_dir, seq_record, options):
    """Set paths, check existing files and run FIMO in parallel on each predicted motif"""
    if not os.path.exists(fimo_dir):
        os.makedirs(fimo_dir)

    args = []
    for plus_minus in os.listdir(meme_dir):
        motif_file = os.path.join(meme_dir, plus_minus, "meme.html")
        sites_file = os.path.join(meme_dir, plus_minus, "binding_sites.fasta")
        output_file = os.path.join(fimo_dir, plus_minus, "fimo.txt")
        output_dir = os.path.join(fimo_dir, plus_minus)

        # input file and binding sites file present and size not zero --> should be fine to use
        if (os.path.isfile(motif_file)
                and os.path.getsize(motif_file) > 0
                and os.path.isfile(sites_file)
                and os.path.getsize(sites_file) > 0):
            # output file already present and size not zero --> do not run FIMO again on this one
            if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
                pass
            else:
                args.append([
                    "fimo",
                    "-verbosity", "1",
                    "-motif", "1",
                    "-thresh", "0.00006",
                    "-oc", output_dir,
                    motif_file,
                    os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta"),
                ])

    errors = run_p(options, args)
    return sum(errors)
