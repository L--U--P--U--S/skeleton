# vim: set fileencoding=utf-8 :
#
# copyright stuff here
# ...
#
#
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# FIXME for sure there is a more clever and elegant way to do this --> don't need this file at all?
"""Run MEME tool, possibly in parallel"""
import sys
import os
from multiprocessing import Pool
import subprocess
# import logging

meme_dir = sys.argv[1]
cpus = int(sys.argv[2])
promoter_sequences = []

def run_meme(fasta):
    # TODO add logging.debug message for each meme call
    exit_code = subprocess.call([
        "meme",
        fasta,
        "-oc", os.path.dirname(fasta) ,
        "-dna",
        "-nostatus",
        "-mod",
        "anr",
        "-nmotifs", "1",
        "-minw", "6",
        "-maxw", "12",
        "-revcomp",
        "-evt", "1.0e+005"])
    return exit_code

for plus_minus in os.listdir(meme_dir):
    input_file = os.path.join(meme_dir, plus_minus, "promoters.fasta")
    output_file = os.path.join(meme_dir, plus_minus, "meme.xml")

    # input file present and size not zero --> should be fine to use
    if os.path.isfile(input_file) and os.path.getsize(input_file) > 0:
        # output file already present and size not zero --> do not run MEME again on this one
        if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
            pass
        else:
            promoter_sequences.append(input_file)

if __name__ == "__main__": # THIS: if clause not "supported" inside the skeleton environment (__init__.py)
    p = Pool(cpus)
    errors = p.map(run_meme, promoter_sequences)

    if sum(errors) > 0: # at least one exit code > 0
        sys.exit(1)
