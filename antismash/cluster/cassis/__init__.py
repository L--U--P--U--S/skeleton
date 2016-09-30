# vim: set fileencoding=utf-8 :
#
# copyright stuff here
# ...
#
#
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Implementation of the CASSIS method for the motif-based prediction of SM gene clusters"""

import logging
import os
import subprocess
import csv
from xml.etree import cElementTree as ElementTree
from pprint import pprint

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from antismash import (
    config,
    utils,
)


name = "cassis"
short_description = "{}: Detect secondary metabolite gene cluster (motif based)".format(name)
priority = 2 # first run hmmdetect plugin to detect core genes (anchor genes) --> seed for cluster prediction with cassis

_required_binaries = [
    ("meme", "4.11.1"),
    ("fimo", "4.11.1"),
]

# all possible promoter sets for motif detection
# plus --> include <plus> promoters downstream the anchor gene's promoter
# minus --> include <minus> promoters upstream the anchor gene's promoter
_plus_minus = []
for plus in range(16):
    for minus in range(16):
        if plus + minus + 1 >= 4:
            _plus_minus.append({"plus": plus, "minus": minus})


def check_prereqs(options):
    """Check for prerequisites"""
    failure_messages = []
    for binary_name, binary_version in _required_binaries:
        if utils.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate executable for {!r}".format(binary_name))
        # TODO: Check binary version here

    return failure_messages


def get_anchor_genes(seq_record):
    """Return all genes which are putative cluster anchor genes"""
    anchor_genes = []

    for feature in seq_record.features:
        if "sec_met" in feature.qualifiers:
            anchor_genes.append(feature.qualifiers["locus_tag"][0])

    return anchor_genes


def get_promoter_id(promoter):
    """Return promoter ID string dependend on involved gene(s)"""
    if len(promoter["id"]) == 1: # 1 gene --> 1 promoter
        return promoter["id"][0]
    else: # 2 bidirectional genes --> 1 shared promoter
        return "{}+{}".format(promoter["id"][0], promoter["id"][1])


def ignore_overlapping(genes):
    """Ignore genes with overlapping locations (skip the second gene of an overlapping couple)"""
    ignored = []

    overlap = True
    while overlap: # check again until we didn't find any overlap in the entire (remaining) gene list
        overlap = False
        non_overlapping = [genes[0]]

        for i in xrange(1, len(genes)):
            # TODO seq_record/genes is sorted by start coordinates? if no: have to consider more cases!
            if (genes[i-1].location.end >= genes[i].location.start) or (genes[i-1].location.start <= genes[i].location.start and genes[i-1].location.end >= genes[i].location.end):
                # A <----->                                             A <---------->
                # B    <----->                                          B   <----->
                logging.warning("Ignoring {} (overlapping with {})".format(utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i-1]))) # TODO info or warning?
                ignored.append(genes[i])
                overlap = True
            else:
                non_overlapping.append(genes[i])

        genes = non_overlapping

    if ignored:
        logging.info("Ignoring {} genes due to overlapping locations".format(len(ignored)))

    return (genes, ignored)


def get_promoters(seq_record, genes, upstream_tss, downstream_tss, options):
    """Compute promoter sequences for each gene in the sequence record"""
    logging.info("Computing promoter sequences")

    min_promoter_length = 6
    max_promoter_length = (upstream_tss + downstream_tss) * 2 + 1

    record_seq_length = len(seq_record.seq)
    promoters = [] # TODO use SeqRecords instead of dicts or even create a class promoter (class cluster, class motif)?
    invalid = 0

    # TODO "format()" or "%s" or "s + s"?
    pos_handle = open(os.path.join(options.outputfoldername, seq_record.name + "_promoter_positions.csv"), "w")
    pos_handle.write("\t".join(["#", "promoter", "start", "end", "length"]) + "\n")
    seq_handle = open(os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta"), "w")

    skip = 0 # helper var for shared promoter of bidirectional genes
    for i in xrange(len(genes)):

        if skip: # two genes share the same promotor --> did computation with first gene, skip second gene
            skip = 0

        elif len(genes) == 1: # only one gene within record

            if genes[i].location.strand == 1:
                if genes[i].location.start - upstream_tss >= 0 and genes[i].location.end > genes[i].location.start + downstream_tss: #1
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i].location.start - upstream_tss < 0 and genes[i].location.end > genes[i].location.start + downstream_tss: #2
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 1, # TODO 0? --> "Note that the start and end location numbering follow Python's scheme, thus a GenBank entry of 123..150 (one based counting) becomes a location of [122:150] (zero based counting)."
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i].location.start - upstream_tss >= 0 and genes[i].location.start + downstream_tss >= genes[i].location.end: #3
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                elif genes[i].location.start - upstream_tss < 0 and genes[i].location.start + downstream_tss >= genes[i].location.end: #7
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 1,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i])) # TODO logging.error --> raise?

            elif genes[i].location.strand == -1:
                if genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= record_seq_length: #4
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss > record_seq_length: #5
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": record_seq_length
                    })
                elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= record_seq_length: #6
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i+1].location.start >= genes[i].location.end - upstream_tss and genes[i].location.end + upstream_tss > record_seq_length: #8
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": record_seq_length
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

        # first gene of the record AND NOT special case #9
        elif ( i == 0 and i != len(genes) ) and not ( genes[i].location.strand == -1 and genes[i+1].location.strand == 1 and genes[i].location.end + upstream_tss >= genes[i+1].location.start - upstream_tss ):

            if genes[i].location.strand == 1:
                if genes[i].location.start - upstream_tss >= 0 and genes[i].location.end > genes[i].location.start + downstream_tss: #1
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i].location.start - upstream_tss < 0 and genes[i].location.end > genes[i].location.start + downstream_tss: #2
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 1,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i].location.start - upstream_tss >= 0 and genes[i].location.start + downstream_tss >= genes[i].location.end: #3
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                elif genes[i].location.start - upstream_tss < 0 and genes[i].location.start + downstream_tss >= genes[i].location.end: #7
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 1,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

            elif genes[i].location.strand == -1:
                if genes[i].location.start < genes[i].location.end - downstream_tss and genes[i+1].location.start > genes[i].location.end + upstream_tss: #4
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i+1].location.start <= genes[i].location.end + upstream_tss: #5
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i+1].location.start - 1
                    })
                elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i+1].location.start > genes[i].location.end + upstream_tss: #6
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i+1].location.start <= genes[i].location.end + upstream_tss and genes[i].location.start >= genes[i].location.end - downstream_tss: #8
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i+1].location.start - 1
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

        # last gene of record
        elif i == len(genes) and not skip:

            if genes[i].location.strand == 1:
                if genes[i-1].location.end < genes[i].location.start - upstream_tss and genes[i].location.end > genes[i].location.start + downstream_tss: #1
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i-1].location.end >= genes[i].location.start - upstream_tss and genes[i].location.end > genes[i].location.start + downstream_tss: #2
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i-1].location.end < genes[i].location.start - upstream_tss and genes[i].location.start + downstream_tss >= genes[i].location.end: #3
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                elif genes[i-1].location.end >= genes[i].location.start - upstream_tss and genes[i].location.start + downstream_tss >= genes[i].location.end: #7
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

            elif genes[i].location.strand == -1:
                if genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= record_seq_length: #4
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss > record_seq_length: #5
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": record_seq_length
                    })
                elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= record_seq_length : #6
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i+1].location.start <= genes[i].location.end + upstream_tss and genes[i].location.end + upstream_tss > record_seq_length: #8
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": record_seq_length
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

        # special-case 9
        elif genes[i].location.strand == -1 and genes[i+1].location.strand == 1 and genes[i].location.end + upstream_tss >= genes[i+1].location.start - upstream_tss:
            if genes[i].location.end > genes[i].location.start + downstream_tss and genes[i].location.start < genes[i].location.end - downstream_tss: #9 (1+4)
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.end - downstream_tss,
                    "end": genes[i+1].location.start + downstream_tss
                })
            elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i+1].location.start + downstream_tss >= genes[i+1].end: #9 (3+4)
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.end - downstream_tss,
                    "end": genes[i+1].end
                })
            elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i+1].end > genes[i+1].location.start + downstream_tss: #9 (1+6)
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.start,
                    "end": genes[i+1].location.start + downstream_tss
                })
            elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i+1].location.start + downstream_tss >= genes[i+1].end: #9 (3+6)
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.start,
                    "end": genes[i+1].end
                })
            else:
                logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

            skip = 1

        # "normal" cases
        elif not skip:

            if genes[i].location.strand == 1:
                if genes[i-1].location.end < genes[i].location.start - upstream_tss and genes[i].location.end > genes[i].location.start + downstream_tss: #1
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i-1].location.end >= genes[i].location.start - upstream_tss and genes[i].location.end > genes[i].location.start + downstream_tss: #2
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.start + downstream_tss
                    })
                elif genes[i-1].location.end < genes[i].location.start - upstream_tss and genes[i].location.start + downstream_tss >= genes[i].location.end: #3
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                elif genes[i-1].location.end >= genes[i].location.start - upstream_tss and genes[i].location.start + downstream_tss >= genes[i].location.end: #7
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

            elif genes[i].location.strand == -1:
                if genes[i].location.start < genes[i].location.end - downstream_tss and genes[i+1].location.start > genes[i].location.end + upstream_tss: #4
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i+1].location.start <= genes[i].location.end + upstream_tss: #5
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i+1].location.start - 1
                    })
                elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i+1].location.start > genes[i].location.end + upstream_tss: #6
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i+1].location.start <= genes[i].location.end + upstream_tss and genes[i].location.start >= genes[i].location.end - downstream_tss: #8
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i+1].location.start - 1
                    })
                else:
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i]))

        # negative start position or stop position "beyond" record --> might happen in very small records
        if promoters[-1]["start"] < 1:
            promoters[-1]["start"] = 1
        if promoters[-1]["end"] > record_seq_length:
            promoters[-1]["end"] = record_seq_length

        # write promoter positions and sequences to file
        if not skip:
            promoter_sequence = seq_record.seq[promoters[-1]["start"]:promoters[-1]["end"]+1]
            promoter_length = len(promoter_sequence)

            invalid_promoter_sequence = ""

            # check if promoter length is valid
            if promoter_length < min_promoter_length or promoter_length > max_promoter_length:
                invalid_promoter_sequence = "length"

            # check if a, c, g and t occur at least once in the promoter sequence
            elif "A" not in promoter_sequence:
                invalid_promoter_sequence = "A"
            elif "C" not in promoter_sequence:
                invalid_promoter_sequence = "C"
            elif "G" not in promoter_sequence:
                invalid_promoter_sequence = "G"
            elif "T" not in promoter_sequence:
                invalid_promoter_sequence = "T"
                # TODO always capital letters?

            if invalid_promoter_sequence:
                invalid += 1

                if invalid_promoter_sequence == "length":
                    logging.warning("Promoter %s is invalid (length is %s)", get_promoter_id(promoters[-1]), promoter_length)
                else:
                    logging.warning("Promoter %s is invalid (sequence without %s)", get_promoter_id(promoters[-1]), invalid_promoter_sequence) # especially SiTaR doesn't like such missings

                # more details for debug logging
                logging.debug("Invalid promoter {}\n start {}\n end {}\n length {}\n".format(promoters[-1]["id"], promoters[-1]["start"], promoters[-1]["end"], promoter_length))

                promoters.pop() # remove last (invalid!) promoter

            else:
                promoters[-1]["seq"] = promoter_sequence

                # write promoter positions to file
                pos_handle.write("\t".join(map(str, [len(promoters), get_promoter_id(promoters[-1]), promoters[-1]["start"], promoters[-1]["end"], promoter_length])) + "\n")

                # write promoter sequences to file
                SeqIO.write(
                    SeqRecord(promoter_sequence, id = get_promoter_id(promoters[-1]),
                    description = "length={}bp".format(promoter_length)),
                    seq_handle,
                    "fasta"
                )

        # check if promoter IDs are unique
        if len(promoters) >= 2 and get_promoter_id(promoters[-1]) == get_promoter_id(promoters[-2]):
            logging.error("Promoter %s occurs at least twice. This may be caused by overlapping gene annotations", get_promoter_id(promoters[-1]))
            # TODO die? raise exception?

    if invalid:
        logging.info("Ignoring {} promoters due to invalid promoter sequences".format(invalid))

    logging.info("Found {} promoter sequences for {} genes".format(len(promoters), len(genes)))
    return promoters


def predict_motifs(anchor, anchor_promoter, promoters, options):
    """Run MEME tool to predict motifs (putative transcription factor binding sites) in promoter sequences"""
    # TODO options --> meme settings
    meme_dir = os.path.join(options.outputfoldername, "meme", anchor) # TODO are anchor gene names save to use for directories?
    motifs = []

    if not os.path.exists(meme_dir):
        os.makedirs(meme_dir)

    # prepare sets of promoter sequences (MEME input)
    indices = set() # helper variable to monitor unique start_index/end_index tupels
    for pm in _plus_minus:
        start_index = anchor_promoter - pm["minus"]
        end_index = anchor_promoter + pm["plus"]

        if start_index < 0: # anchor promoter near beginning of record --> truncate
            logging.debug("Promoter set +{}_-{} exceeds upstream record border".format(pm["plus"], pm["minus"]))
            start_index = 0

        if end_index > len(promoters) - 1: # anchor promoter near end of record --> truncate
            logging.debug("Promoter set +{}_-{} exceeds downstream record border".format(pm["plus"], pm["minus"]))
            end_index = len(promoters) - 1

        # discard promoter sets, which reappear due to truncation
        if (start_index, end_index) not in indices:
            indices.add((start_index, end_index))
            motifs.append({"plus": pm["plus"], "minus": pm["minus"], "score": None})

            pm_dir = os.path.join(meme_dir, "+{}_-{}".format(pm["plus"], pm["minus"]))
            if not os.path.exists(pm_dir):
                os.makedirs(pm_dir)

            # write promoter sequences to fasta file, in respective "plus-minus" subdir
            with open(os.path.join(pm_dir, "promoters.fasta"), "w") as pm_handle:
                for i in xrange(start_index, end_index + 1):
                    seq = SeqRecord(promoters[i]["seq"], id = get_promoter_id(promoters[i]), description = "length={}bp".format(len(promoters[i]["seq"])))
                    if i == anchor_promoter: # mark anchor gene (must be part of id, otherwise MEME woun't recognize it)
                        # seq.description += " ANCHOR"
                        seq.id += "__ANCHOR"
                    SeqIO.write(seq, pm_handle, "fasta")
        else:
            logging.debug("Duplicate promoter set +{}_-{}".format(pm["plus"], pm["minus"]))

    # run MEME
    # FIXME for sure there is a more clever and elegant way to do this
    exit_code = subprocess.call(["python", os.path.join(os.path.dirname(os.path.realpath(__file__)), "meme.py"), meme_dir, str(options.cpus)])
    if exit_code != 0:
        logging.error("meme.py discovered a problem (exit code {})".format(exit_code))

    # analyse MEME results
    for motif in motifs:
        xml_file = os.path.join(meme_dir, "+{}_-{}".format(motif["plus"], motif["minus"]), "meme.xml")
        e = ElementTree.parse(xml_file).getroot()
        reason = e.find("model/reason_for_stopping").text
        anchor_seq_id = ""

        # no motif found for given e-value cutoff :-(
        if "Stopped because motif E-value > " in reason:
            logging.debug("Motif +{:02d}/-{:02d}: e-value exceeds cutoff".format(motif["plus"], motif["minus"]))

        # motif(s) found :-)
        elif "Stopped because requested number of motifs (1) found" in reason:

            # find anchor genes' sequence_id
            training_set = e.findall("training_set/sequence") # all promoter sequences passed to MEME
            for i in xrange(len(training_set)):
                if "__ANCHOR" in training_set[i].attrib["name"]:
                    anchor_seq_id = training_set[i].attrib["id"] # e.g. id=sequence_1

            # only accept motifs which occur in the anchor genes promoter
            contributing_sites = e.findall("motifs/motif/contributing_sites/contributing_site") # sequences which contributed to the motif
            if anchor_seq_id in map(lambda site: site.attrib["sequence_id"], contributing_sites): # could do it with filter(), but this seems to be less clear
                # save motif score
                motif["score"] = e.find("motifs/motif").attrib["e_value"] # one motif, didn't ask MEME for more

                # save sequence sites which represent the motif
                motif["seqs"] = ["".join(map(lambda letter: letter.attrib["letter_id"], site.findall("site/letter_ref"))) for site in contributing_sites]

                # write sites to fasta file
                with open(os.path.join(meme_dir, "+{}_-{}".format(motif["plus"], motif["minus"]), "binding_sites.fasta"), "w") as handle:
                    handle.write(">{}__+{}_-{}\n".format(anchor, motif["plus"], motif["minus"]))
                    handle.write("\n".join(motif["seqs"]))

                logging.debug("Motif +{:02d}/-{:02d}: e-value = {}".format(motif["plus"], motif["minus"], motif["score"]))
            else:
                logging.debug("Motif +{:02d}/-{:02d}: does not occur in anchor gene promoter".format(motif["plus"], motif["minus"]))

        # unexpected reason, don't know why MEME stopped :-$
        else:
            logging.error("MEME stopped unexpectedly (reason: {})".format(reason))

    return filter(lambda m: m["score"] is not None, motifs)


def search_motifs(anchor, anchor_promoter, motifs, promoters, seq_record, options):
    """Run FIMO tool to find occurrences of previously predicted motifs (by MEME)"""
    meme_dir = os.path.join(options.outputfoldername, "meme", anchor)
    fimo_dir = os.path.join(options.outputfoldername, "fimo", anchor) # TODO are anchor gene names save to use for directories?
    promoter_sequences = os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta")

    if not os.path.exists(fimo_dir):
        os.makedirs(fimo_dir)

    # run FIMO
    # FIXME for sure there is a more clever and elegant way to do this
    exit_code = subprocess.call(["python", os.path.join(os.path.dirname(os.path.realpath(__file__)), "fimo.py"), meme_dir, fimo_dir, promoter_sequences, str(options.cpus)])
    if exit_code != 0:
        logging.error("fimo.py discovered a problem (exit code {})".format(exit_code))

    # analyse FIMO results
    for motif in motifs:
        motif["hits"] = {}
        with open(os.path.join(fimo_dir, "+{}_-{}".format(motif["plus"], motif["minus"]), "fimo.txt"), "r") as handle: # TODO r or rb?
            table = csv.reader(handle, delimiter = "\t")
            for row in table:
                if not row[0].startswith("#"): # skip comment lines
                    seq_id = row[1]
                    if seq_id in motif["hits"]:
                        motif["hits"][seq_id] += 1
                    else:
                        motif["hits"][seq_id] = 1

        # write binding sites per promoter to file
        with open(os.path.join(fimo_dir, "+{}_-{}".format(motif["plus"], motif["minus"]), "bs_per_promoter.csv"), "w") as handle: # TODO w or wb?
            table = csv.writer(handle, delimiter = "\t", lineterminator = "\n")
            table.writerow(["#", "promoter", "binding sites"]) # table head
            for i in xrange(len(promoters)):
                promoter = get_promoter_id(promoters[i])
                if promoter in motif["hits"]:
                    table.writerow([i+1, promoter, motif["hits"][promoter]])
                else:
                    table.writerow([i+1, promoter, 0])

        percentage = float(len(motif["hits"])) / float(len(promoters)) * 100 # float!
        if percentage == 0.0:
            # too low
            logging.debug("Motif +{:02d}/-{:02d}: occurs in {} promoters (no hits)".format(motif["plus"], motif["minus"], len(motif["hits"])))
            motif["hits"] = None
        elif percentage > 14.0: # TODO options
            # too high
            logging.debug("Motif +{:02d}/-{:02d}: occurs in {} promoters, {:.2f}% of all promoters (too many)".format(motif["plus"], motif["minus"], len(motif["hits"]), percentage))
            motif["hits"] = None
        elif get_promoter_id(promoters[anchor_promoter]) not in motif["hits"]: # not in achor promoter
            # no site in anchor promoter
            logging.debug("Motif +{:02d}/-{:02d}: not hits in the promoter of the anchor gene".format(motif["plus"], motif["minus"]))
            motif["hits"] = None
        else:
            # everything ok
            logging.debug("Motif +{:02d}/-{:02d}: occurs in {} promoters, {:.2f}% of all promoters".format(motif["plus"], motif["minus"], len(motif["hits"]), percentage))

    return filter(lambda m: m["hits"] is not None, motifs)


def find_islands(anchor_promoter, motifs, promoters, options):
    """Find islands of binding sites (previously found by FIMO) around anchor gene to define cluster borders"""
    max_gap_length = 2 # TODO options
    islands = []

    for motif in motifs:
        # create list with binding sites per promoter
        bs_per_promoter = [0] * len(promoters) # first set number of binding sites to 0
        for i in xrange(len(promoters)):
            if get_promoter_id(promoters[i]) in motif["hits"]: # second set actual number of binding sites, if any
                bs_per_promoter[i] = motif["hits"][get_promoter_id(promoters[i])]

        # upstream
        start = anchor_promoter # init upstream cluster border
        i = anchor_promoter # init position of anchor gene's promoter
        while i > 0:

            # promoter with binding site
            # … 1 …
            if bs_per_promoter[i-1] >= 1:
                start -= 1

            # no binding site, gap with length 1
            # … 1 0 1  …
            elif (i - 2 >= 0
                  and max_gap_length >= 1
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] >= 1):
                start -= 2
                i -= 1

            # no binding site, gap with length 2
            # … 1 0 0 1  …
            elif (i - 3 >= 0
                  and max_gap_length >= 2
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] >= 1):
                start -= 3
                i     -= 2

            # no binding site, gap with length 3
            # … 1 0 0 0 1  …
            elif (i - 4 >= 0
                  and max_gap_length >= 3
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] == 0
                  and bs_per_promoter[i-4] >= 1):
                start -= 4
                i     -= 3

            # no binding site, gap with length 4
            # … 1 0 0 0 0 1  …
            elif (i - 5 >= 0
                  and max_gap_length >= 4
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] == 0
                  and bs_per_promoter[i-4] == 0
                  and bs_per_promoter[i-5] >= 1):
                start -= 5
                i     -= 4

            # no binding site, gap with length 5
            # … 1 0 0 0 0 0 1  …
            elif (i - 3 >= 6
                  and max_gap_length >= 5
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] == 0
                  and bs_per_promoter[i-4] == 0
                  and bs_per_promoter[i-2] == 5
                  and bs_per_promoter[i-6] >= 1):
                start -= 6
                i     -= 5

            # gap too long, stop upstream cluster extension
            else:
                break

            i -= 1

        # downstream
        i = anchor_promoter # reset position of anchor gene's promoter
        end = anchor_promoter # init downstream cluster border
        while i < len(bs_per_promoter) - 1:

            # promoter with binding site(s)
            # … 1 …
            if bs_per_promoter[i+1] > 0:
                end += 1

            # no binding site, gap with length 1
            # … 1 0 1  …
            elif (i + 2 < len(bs_per_promoter)
                  and max_gap_length >= 1
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] >= 1):
                end += 2
                i += 1

            # no binding site, gap with length 2
            # … 1 0 0 1  …
            elif (i + 3 < len(bs_per_promoter)
                  and max_gap_length >= 2
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] >= 1):
                end += 3
                i += 2

            # no binding site, gap with length 3
            # … 1 0 0 0 1  …
            elif (i + 4 < len(bs_per_promoter)
                  and max_gap_length >= 3
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] == 0
                  and bs_per_promoter[i+4] >= 1):
                end += 4
                i += 3

            # no binding site, gap with length 4
            # … 1 0 0 0 0 1  …
            elif (i + 5 < len(bs_per_promoter)
                  and max_gap_length >= 4
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] == 0
                  and bs_per_promoter[i+4] == 0
                  and bs_per_promoter[i+5] >= 1):
                end += 5
                i += 4

            # no binding site, gap with length 5
            # … 1 0 0 0 0 0 1  …
            elif (i + 6 < len(bs_per_promoter)
                  and max_gap_length >= 5
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] == 0
                  and bs_per_promoter[i+4] == 0
                  and bs_per_promoter[i+5] == 0
                  and bs_per_promoter[i+6] >= 1):
                end += 6
                i += 5

            # gap too long, stop downstream cluster extension
            else:
                break

            i += 1

        logging.debug("Motif +{:02d}/-{:02d}: island [{}, {}]".format(motif["plus"], motif["minus"], get_promoter_id(promoters[start]), get_promoter_id(promoters[end])))
        islands.append({"start": promoters[start], "end": promoters[end], "motif": motif}) # TODO save promoters or promoter ids or first/last gene?

    return islands


def sort_by_abundance(islands):
    """sort upstream (start) and downstream (end) borders of islands by abundance"""
    # count border abundance
    # start/end are treated independently!
    starts = {}
    ends = {}
    for i in islands:
        if i["start"]["id"][0] in starts:
            starts[i["start"]["id"][0]]["abund"] += 1
        else:
            starts[i["start"]["id"][0]] = {"abund": 1}

        if i["end"]["id"][-1] in ends:
            ends[i["end"]["id"][-1]]["abund"] += 1
        else:
            ends[i["end"]["id"][-1]] = {"abund": 1}

        # also keep track of motif score
        # --> sort by score if same abundance occurs more than once
        starts[i["start"]["id"][0]]["mscore"] = i["motif"]["score"]
        ends[i["end"]["id"][-1]]["mscore"] = i["motif"]["score"]

    # compute sum of start and end abundance, remove duplicates, sort descending
    abundances_sum_sorted = sorted(set([s["abund"] + e["abund"] for s in starts.values() for e in ends.values()]), reverse=True)
    # compute sum of start and end motif score, remove duplicates, sort ascending
    scores_sum_sorted = sorted(set([float(s["mscore"]) + float(e["mscore"]) for s in starts.values() for e in ends.values()]))
    # sort by value (=abundance) of start, descending
    starts_sorted = sorted(starts, key=starts.get, reverse=True)
    # sort by value (=abundance) of end, descending
    ends_sorted = sorted(ends, key=ends.get, reverse=True)

    clusters = []
    for abundance in abundances_sum_sorted:
        # list from highest (best) to lowest (worst) abundance
        for score in scores_sum_sorted:
            # list from lowest (best) to highest (worst) motif score/e-value
            for start in starts_sorted:
                for end in ends_sorted:
                    if starts[start]["abund"] + ends[end]["abund"] == abundance and float(starts[start]["mscore"]) + float(ends[end]["mscore"]) == score:
                        logging.debug("Abundance {}, score {:.1e}: [{}, {}, {} -- {}, {}, {}]".format(abundance, score, start, starts[start]["abund"], starts[start]["mscore"], end, ends[end]["abund"], ends[end]["mscore"]))
                        clusters.append( {"start": {"gene": start, "abundance": starts[start]["abund"], "score": starts[start]["mscore"]}, "end": {"gene": end, "abundance": ends[end]["abund"], "score": ends[end]["mscore"]}})
                        # list[ dict{left border: dict, right border: dict} ]

    return clusters



def detect(seq_record, options):
    """Use core genes (anchor genes) from hmmdetect as seeds to detect gene clusters"""
    logging.info("Detecting gene clusters using CASSIS method")

    # TODO options? cassis settings/parameters?

    # get core genes from hmmdetect --> necessary CASSIS input, aka "anchor genes"
    anchor_genes = get_anchor_genes(seq_record)

    ### TESTING ###
    anchor_genes.append("AFUA_6G00240")
    ### TESTING ###

    # compute promoter sequences/regions --> necessary for motif prediction (MEME and FIMO input)
    upstream_tss = 1000; # nucleotides upstream TSS
    downstream_tss = 50; # nucleotides downstream TSS
    genes, ignored_genes = ignore_overlapping(utils.get_all_features_of_type(seq_record, "gene"))
    promoters = get_promoters(seq_record, genes, upstream_tss, downstream_tss, options)

    if len(promoters) < 3:
        logging.warning("Sequence {!r} yields less than 3 promoter regions, skipping cluster detection".format(seq_record.name))
    else:
        if len(promoters) < 40:
            logging.warning("Sequence {!r} yields only {} promoter regions. Cluster detection on small sequences may lead to incomplete cluster predictions".format(seq_record.name, len(promoters)))

        for anchor in anchor_genes:
            logging.info("Detecting cluster around anchor gene {}".format(anchor))

            anchor_promoter = None
            for i in xrange(len(promoters)):
                # the promoter ID string is not equal to the anchor ID string! (it's somewhere "in between")
                if anchor in promoters[i]["id"]:
                    anchor_promoter = i
                    break

            if anchor_promoter is None:
                logging.warning("No promoter region for {}, skipping anchor gene".format(anchor))
                continue

            # predict motifs with MEME around the anchor gene
            motifs = predict_motifs(anchor, anchor_promoter, promoters, options)

            if len(motifs) == 0:
                logging.info("Could not predict motifs around {}, skipping anchor gene".format(anchor))
                continue

            # search predicted binding sites with FIMO in all promoter sequences
            # and count number of occurrences per promoter
            motifs = search_motifs(anchor, anchor_promoter, motifs, promoters, seq_record, options)

            if len(motifs) == 0:
                logging.info("Could not find motif occurrences for {}, skipping anchor gene".format(anchor))
                continue

            # find islands of binding sites around anchor gene
            islands = find_islands(anchor_promoter, motifs, promoters, options)
            logging.debug("{} cluster predictions for {}".format(len(islands), anchor))

            # return cluster predictions sorted by border abundance
            # most abundant --> "best" prediction --> index 0
            # so far, only use best prediction (this may change in later versions of this plugin)
            cluster_predictions = sort_by_abundance(islands)
            logging.info("Best prediction: {} -- {}".format(cluster_predictions[0]["start"]["gene"], cluster_predictions[0]["end"]["gene"]))

            # warn if cluster prediction right at or next to record (~ contig) border
            all_genes = map(lambda g: utils.get_gene_id(g), utils.get_all_features_of_type(seq_record, "gene"))
            start_index = all_genes.index(cluster_predictions[0]["start"]["gene"])
            end_index = all_genes.index(cluster_predictions[0]["end"]["gene"])
            if start_index < 10:
                logging.warning("Upstream cluster border located at or next to sequence record border, prediction could have been truncated by record border")
            if end_index > len(all_genes) - 10:
                logging.warning("Downstream cluster border located at or next to sequence record border, prediction could have been truncated by record border")

            # warn if cluster prediction too short (includes less than 3 genes)
            if end_index - start_index + 1 < 3:
                logging.warning("Cluster is very short (less than 3 genes). Prediction may be questionable.")

            # warn if ignored gene (overlapping with anthor gene, see ignore_overlapping()) would have been part of the cluster
            for ignored_gene in map(lambda g: utils.get_gene_id(g), ignored_genes):
                if ignored_gene in all_genes[start_index : end_index + 1]:
                    logging.warning("Ignored gene {} could have affected the prediction".format(ignored_gene))
                    break


def get_versions(options):
    """Get all utility versions"""
    # TODO meme, fimo, …
    return []


# TODO implement for cassis?
#def store_detection_details(rulesdict, seq_record):
#    '''Store the details about why a cluster was detected'''
#    clusters = utils.get_cluster_features(seq_record)
#    for cluster in clusters:
#        type_combo = utils.get_cluster_type(cluster)
#        if '-' in type_combo:
#            clustertypes = type_combo.split('-')
#        else:
#            clustertypes = [type_combo]
#
#        if not 'note' in cluster.qualifiers:
#            cluster.qualifiers['note'] = []
#        rule_string = "Detection rule(s) for this cluster type:"
#        for clustertype in clustertypes:
#            rule_string += " %s: (%s);" % (clustertype, rulesdict[clustertype][0])
#
#        cluster.qualifiers['note'].append(rule_string)

# TODO implement for cassis? (Kai?)
#def _update_sec_met_entry(feature, results, clustertype, nseqdict):
#    '''Update the sec_met entry for a feature'''
#    result = "; ".join(["%s (E-value: %s, bitscore: %s, seeds: %s)" % (
#        res.query_id, res.evalue, res.bitscore, nseqdict.get(res.query_id, '?')) for res in results])
#
#    if not 'sec_met' in feature.qualifiers:
#        feature.qualifiers['sec_met'] = [
#            "Type: %s" % clustertype,
#            "Domains detected: %s" % (result),
#            "Kind: biosynthetic"
#        ]
#    else:
#        for ann in feature.qualifiers['sec_met']:
#            if not ann.startswith("Domains detected"):
#                continue
#            ann += "Domains detected: %s" % (result)
