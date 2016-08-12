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
"""Run FIMO tool, possibly in parallel"""
import sys
import os
from multiprocessing import Pool
import subprocess

meme_dir = sys.argv[1]
fimo_dir = sys.argv[2]
promoter_sequences = sys.argv[3]
cpus = int(sys.argv[4])
motifs = []

def run_fimo(motif):
    exit_code = subprocess.call(["fimo", "-verbosity", "1", "-motif", "1", "-thresh", "0.00006", "-oc", motif["out_dir"], motif["motif"], promoter_sequences])
    return exit_code

for plus_minus in os.listdir(meme_dir):
    motif_file = os.path.join(meme_dir, plus_minus, "meme.html")
    sites_file = os.path.join(meme_dir, plus_minus, "binding_sites.fasta")
    output_file = os.path.join(fimo_dir, plus_minus, "fimo.txt")
    output_dir = os.path.join(fimo_dir, plus_minus)

    # input file and binding sites file present and size not zero --> should be fine to use
    if os.path.isfile(motif_file) and os.path.getsize(motif_file) > 0 and os.path.isfile(sites_file) and os.path.getsize(sites_file) > 0:
        # output file already present and size not zero --> do not run FIMO again on this one
        if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
            pass
        else:
            motifs.append({"motif": motif_file, "out_dir": output_dir})

if __name__ == "__main__": # THIS: if clause not "supported" inside the skeleton environment (__init__.py)
    p = Pool(cpus)
    errors = p.map(run_fimo, motifs)

    if sum(errors) > 0: # at least one exit code > 0
        sys.exit(1)
