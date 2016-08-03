# vim: set fileencoding=utf-8 :
#
# copyright stuff here
# ...
#
#
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Load sequences from a local file"""
import logging
import Bio
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
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


def check_prereqs(options):
    """Check for prerequisites"""
    failure_messages = []
    for binary_name, binary_version in _required_binaries:
        if utils.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate executable for {!r}".format(binary_name))
        # TODO: Check binary version here

    return failure_messages


def get_all_anchor_genes(seq_record):
    """Return all genes which are putative cluster anchor genes"""
    anchor_genes = []

    for feature in seq_record.features:
        if "sec_met" in feature.qualifiers:
            anchor_genes.append(feature.qualifiers["locus_tag"][0])

    return anchor_genes


def get_promoter_id(promoter):
    if len(promoter["id"]) == 1: # 1 gene --> 1 promoter
        return promoter["id"][0]
    else: # 2 bidirectional genes --> 1 shared promoter
        return "{}+{}".format(promoter["id"][0], promoter["id"][1])


def ignore_overlapping_genes(genes):
    """Ignore genes with overlapping locations (skip the second one of an overlapping couple)"""

    ignored = 0

    overlap = True
    while overlap: # check again until we didn't find any overlap in the entire (remaining) gene list
        overlap = False
        non_overlapping = [genes[0]]

        for i in xrange(1, len(genes)):
            # TODO seq_record/genes is sorted by start coordinates? if no: have to consider more cases!
            if (genes[i-1].location.end >= genes[i].location.start) or (genes[i-1].location.start <= genes[i].location.start and genes[i-1].location.end >= genes[i].location.end):
                # A <----->                                             A <---------->
                # B    <----->                                          B   <----->
                logging.warning("Ignoring {} (overlapping with {})".format(utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i-1]))) # # TODO info or warning?
                ignored += 1
                overlap = True
            else:
                non_overlapping.append(genes[i])

        genes = non_overlapping

    if ignored:
        logging.info("Ignoring {} genes due to overlapping locations".format(ignored))

    return genes


def get_all_promoters(seq_record, upstream_tss, downstream_tss, options):
    """Compute promoter sequences for each gene"""
    logging.info("Computing promoter sequences")

    min_promoter_length = 6
    max_promoter_length = (upstream_tss + downstream_tss) * 2 + 1

    genes = ignore_overlapping_genes(utils.get_all_features_of_type(seq_record, "gene"))
    contig_length = len(seq_record.seq) # TODO is 1 seq_record == 1 contig always true?
    promoters = []
    invalid = 0

    # TODO "format()" or "%s" or "s + s"?
    pos_handle = open("{}/{}_promoter_positions.csv".format(options.outputfoldername, seq_record.name), "w")
    pos_handle.write("\t".join(["#", "promoter", "start", "end", "length"]) + "\n")
    seq_handle = open("{}/{}_promoter_sequences.fasta".format(options.outputfoldername, seq_record.name), "w")

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
                    logging.error("Problem with promoter of gene '%s'", utils.get_gene_id(genes[i])) # TODO logging.error --> die?

            elif genes[i].location.strand == -1:
                if genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= contig_length: #4
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss > contig_length: #5
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": contig_length
                    })
                elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= contig_length: #6
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i+1].location.start >= genes[i].location.end - upstream_tss and genes[i].location.end + upstream_tss > contig_length: #8
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": contig_length
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
                if genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= contig_length: #4
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i].location.start < genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss > contig_length: #5
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": contig_length
                    })
                elif genes[i].location.start >= genes[i].location.end - downstream_tss and genes[i].location.end + upstream_tss <= contig_length : #6
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                elif genes[i+1].location.start <= genes[i].location.end + upstream_tss and genes[i].location.end + upstream_tss > contig_length: #8
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": contig_length
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
        if promoters[-1]["end"] > contig_length:
            promoters[-1]["end"] = contig_length

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

                promoters.pop() # remove last (invalid!) promoter

            else:
                # write promoter positions to file
                pos_handle.write("\t".join(map(str, [len(promoters), get_promoter_id(promoters[-1]), promoters[-1]["start"], promoters[-1]["end"], promoter_length])) + "\n")

                # write promoter sequences to file
                SeqIO.write(SeqRecord(promoter_sequence, id = get_promoter_id(promoters[-1]), description = "length={}bp".format(promoter_length)), seq_handle, "fasta")

        # check if promoter IDs are unique
        if len(promoters) >= 2 and get_promoter_id(promoters[-1]) == get_promoter_id(promoters[-2]):
            logging.error("Promoter %s occurs at least twice. This may be caused by overlapping gene annotations", get_promoter_id(promoters[-1]))
            # TODO die? raise exception?

    if invalid:
        logging.info("Ignoring {} promoters due to invalid promoter sequences".format(invalid))

    logging.info("Found {} promoter sequences for {} genes".format(len(promoters), len(genes)))
    return promoters


def detect(seq_record, options):
    """Predict SM gene clusters using CASSIS (cluster assignment by islands of sites)"""
    logging.info("Detecting gene clusters using CASSIS method")

    # TODO options? cassis settings/parameters?

    # print dir(seq_record)
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__iter__', '__le___', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'upper']

    # get core genes from hmmdetect --> necessary CASSIS input, aka "anchor genes"
    anchor_genes = get_all_anchor_genes(seq_record)

    # compute promoter sequences/regions --> necessary for motif prediction (MEME and FIMO input)
    upstream_tss = 1000; # nucleotides upstream TSS
    downstream_tss = 50; # nucleotides downstream TSS
    promoters = get_all_promoters(seq_record, upstream_tss, downstream_tss, options)

    if len(promoters) < 40:
        logging.warning("Sequence {!r} yields only {} promoter regions. Cluster detection on small sequences may lead to incomplete cluster predictions".format(seq_record.name, len(promoters)))

    if len(promoters) < 3:
        logging.error("Sequence {!r} yields less than 3 promoter regions. Skipping cluster detection".format(seq_record.name))
    else:
        find_clusters(promoters, anchor_genes, options)


def get_versions(options):
    """Get all utility versions"""
    # TODO meme, fimo, …
    return []


def find_clusters(promoters, anchor_genes, options):
    """Use core genes (anchor genes) as seeds to detect gene clusters"""

    for anchor in anchor_genes:
        logging.info("Detecting cluster around anchor gene {}".format(anchor))

        anchor_promoter = None
        for i in xrange(0, len(promoters)):
            # the promoter ID string is not equal to the anchor ID string! (it's somewhere "in between")
            # ----> TODO save gene(s!) of promoter region as list as promoter id? <-------
            if anchor in promoters[i]["id"]:
                anchor_promoter = i

        if not anchor_promoter:
            logging.warning("No promoter region for {}, skipping".format(anchor))




#    #Functions that detects the gene clusters based on the identified core genes
#    features = utils.get_cds_features(seq_record)
#    clustertype = ""
#    clusters = []
#    cfg = config.get_config()
#    clusternr = cfg.next_clusternr
#    state = "seed"
#    last_cluster_extended = False
#
#    for feature in features:
#        if state == "seed":
#            if not 'sec_met' in feature.qualifiers or len([feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat]) == 0:
#                continue
#            coregenetype = [feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat][0].partition("Type: ")[2]
#            clustertype = coregenetype
#            if "-" in clustertype:
#                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
#                extension = max([rulesdict[value][2] for value in clustertype.split("-")])
#            else:
#                cutoff = rulesdict[clustertype][1]
#                extension = rulesdict[clustertype][2]
#            start = max(feature.location.start - extension, 0)
#            if int(start) < int(feature.location.end):
#                loc = FeatureLocation(start, feature.location.end)
#            else:
#                loc = FeatureLocation(feature.location.end, start)
#            new_cluster = SeqFeature(loc, type="cluster")
#            new_cluster.qualifiers['product'] = [clustertype]
#            new_cluster.qualifiers['note'] = ["Cluster number: " + str(clusternr)]
#            clusternr += 1
#            clusters.append(new_cluster)
#            state = "extend"
#            last_cluster_extended = False
#
#        elif state == "extend":
#            cluster = clusters[-1]
#            if "-" in clustertype:
#                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
#                extension = max([rulesdict[value][2] for value in clustertype.split("-")])
#            else:
#                cutoff = rulesdict[clustertype][1]
#                extension = rulesdict[clustertype][2]
#            if not 'sec_met' in feature.qualifiers:
#                if feature.location.start > cluster.location.end + cutoff:
#                    # ok, no hits for too long, done with this cluster
#                    state = "seed"
#                    end = min(len(seq_record), cluster.location.end + extension)
#                    if int(cluster.location.start) < int(end):
#                        cluster.location = FeatureLocation(cluster.location.start, end)
#                    else:
#                        cluster.location = FeatureLocation(end, cluster.location.start)
#                    last_cluster_extended = True
#                continue
#            coregenetype = [feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat][0].partition("Type: ")[2]
#            if coregenetype not in clustertype and coregenetype != 'other':
#                if clustertype == "other" and coregenetype != "other":
#                    clustertype = coregenetype
#                else:
#                    if "-" not in coregenetype:
#                        if coregenetype not in clustertype and coregenetype != "other":
#                            clustertype += "-" + coregenetype
#                    else:
#                        for partctype in [ctype for ctype in coregenetype.split("-") if ctype != "other"]:
#                            if partctype not in clustertype:
#                                clustertype += "-" + partctype
#                cluster.qualifiers['product'] = [clustertype]
#            if int(cluster.location.start) < int(feature.location.end):
#                cluster.location = FeatureLocation(cluster.location.start, feature.location.end)
#            else:
#                cluster.location = FeatureLocation(feature.location.end, cluster.location.start)
#
#    if len(clusters) > 0 and not last_cluster_extended:
#        cluster = clusters[-1]
#        end = min(len(seq_record), feature.location.end + extension)
#        if int(cluster.location.start) < int(end):
#            cluster.location = FeatureLocation(cluster.location.start, end)
#        else:
#            cluster.location = FeatureLocation(end, cluster.location.start)
#
#    seq_record.features.extend(clusters)
#    cfg.next_clusternr = clusternr


#def filter_results(results, results_by_id):
#    #Filter results by comparing scores of different models (for PKS systems)
#    for line in open(utils.get_full_path(__file__, "filterhmmdetails.txt"), "r").read().split("\n"):
#        filterhmms = line.split(",")
#        for cds in results_by_id.keys():
#            cdsresults = results_by_id[cds]
#            hmmhits = [hit.query_id for hit in cdsresults]
#            #Check if multiple competing HMM hits are present
#            competing_hits = set(hmmhits) & set(filterhmms)
#            if len(competing_hits) > 1:
#                #Identify overlapping hits
#                overlapping_groups = []
#                for hit in cdsresults:
#                    for otherhit in [cdsresult for cdsresult in cdsresults if hit != cdsresult]:
#                        overlap = len(set(range(hit.hit_start, hit.hit_end)) & set(range(otherhit.hit_start, otherhit.hit_end)))
#                        if overlap > 20:
#                            added = "n"
#                            for group in overlapping_groups:
#                                if hit in group and otherhit in group:
#                                    added = "y"
#                                    break
#                                elif hit in group and otherhit not in group:
#                                    group.append(otherhit)
#                                    added = "y"
#                                    break
#                                elif hit not in group and otherhit in group:
#                                    group.append(hit)
#                                    added = "y"
#                                    break
#                            if added == "n":
#                                overlapping_groups.append([hit, otherhit])
#                #Remove worst-scoring of overlapping hits
#                for group in overlapping_groups:
#                    highestscore = max([hit.bitscore for hit in group])
#                    hit_with_highestscore = group[[hit.bitscore for hit in group].index(highestscore)]
#                    to_delete = [hit for hit in group if hit != hit_with_highestscore]
#                    for res in [res for res in results]:
#                        if res in to_delete:
#                            del results[results.index(res)]
#                            del results_by_id[cds][results_by_id[cds].index(res)]
#    return results, results_by_id


#def create_rules_dict():
#    "Create a cluster rules dictionary from the cluster rules file"
#    rulesdict = {}
#    first = True
#    #TODO: We should move all user-customizable files into config subdirectory; the rulefiles are redundant also in hmm_detection_dblookup
#    for line in open(utils.get_full_path(__file__, "cluster_rules.txt"), "r"):
#        # skip the first line with the legend
#        if first:
#            first = False
#            continue
#        parts = line.split('\t')
#        if len(parts) < 3:
#            continue
#        key = parts.pop(0)
#        rules = parts.pop(0)
#        cutoff = int(parts.pop(0)) * 1000
#        extension = int(parts.pop(0)) * 1000
#        rulesdict[key] = (rules, cutoff, extension)
#    return rulesdict


#def apply_cluster_rules(results_by_id, feature_by_id, rulesdict, inclusive):
#    "Apply cluster rules to determine if HMMs lead to secondary metabolite core gene detection"
#    typedict = {}
#    for cds in results_by_id.keys():
#        _type = "none"
#        cdsresults = [res.query_id for res in results_by_id[cds]]
#        for clustertype in rulesdict.keys():
#            if not inclusive and clustertype.startswith('cf_'):
#                continue
#            single_rules = [rule for rule in rulesdict[clustertype][0].split(" or ") if " & " not in rule and "cluster(" not in rule]
#            combined_rules = [rule[1:-1] for rule in rulesdict[clustertype][0].split(" or ") if " & " in rule]
#            cluster_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "cluster(" in rule]
#            minimum_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "minimum(" in rule]
#            if "-" in clustertype:
#                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
#            else:
#                cutoff = rulesdict[clustertype][1]
#            #Assign cluster type if a single argument rule matches
#            #Example rule format: "Domain1"
#            if len(set(single_rules) & set(cdsresults)) >= 1:
#                if not (_type != "none" and clustertype == "other"):
#                    if _type == "none" or _type == "other" or _type == clustertype:
#                        _type = clustertype
#                    elif clustertype not in _type:
#                        _type = clustertype + "-" + _type
#                if _type != "other":
#                    continue
#            #Assign cluster type if a combinatorial argument rule matches
#            #Example rule format: "(Domain1 & Domain2)"
#            for rule in combined_rules:
#                required_matches = rule.split(" & ")
#                if len(set(required_matches) & set(cdsresults)) == len(required_matches):
#                    if not (_type != "none" and clustertype == "other"):
#                        if _type == "none" or _type == "other" or _type == clustertype:
#                            _type = clustertype
#                        elif clustertype not in _type:
#                            _type = clustertype + "-" + _type
#            if _type == clustertype and _type != "other":
#                continue
#            #Assign cluster type if distance-based combinatorial parameter matches
#            #Example rule format: "cluster(Domain1,Domain2)"
#            for rule in cluster_rules:
#                #Find which needed domains are already found in present CDS
#                required_matches = rule.split(",")
#                cluster_results = list(set(required_matches) & set(cdsresults))
#                missing_results = list(set(required_matches) - set(cdsresults))
#                #If more than one, search nearby CDSs for the other needed domains
#                if len(cluster_results) > 0:
#                    locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
#                    for othercds in results_by_id.keys():
#                        needed_hits_found = list(set([res.query_id for res in results_by_id[othercds]]) & set(missing_results))
#                        if len(needed_hits_found) > 0:
#                            feature = feature_by_id[othercds]
#                            flocations = [feature.location.start, feature.location.end]
#                            #If hit found in nearby CDS, add it to the set of found domains relevant to the present rule
#                            if min([abs(max(locations) - min(flocations)), abs(max(locations) - min(flocations))]) < cutoff:
#                                cluster_results.extend(needed_hits_found)
#                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
#                if len(list(set(required_matches) - set(cluster_results))) == 0:
#                    if not (_type != "none" and clustertype == "other"):
#                        if _type == "none" or _type == "other" or _type == clustertype:
#                            _type = clustertype
#                        elif clustertype not in _type:
#                            _type = clustertype + "-" + _type
#                        break
#            #Assign cluster type if distance-based combinatorial parameter matches with minimum number of domain hits
#            #Example rule format: "minimum(4,[Domain1,Domain2], [Domain1])"
#            #This denotes that at least 4 hits of either Domain1 or Domain2 should be found in the same region with at least 1 hit from Domain 1
#            for rule in minimum_rules:
#                #Find which needed domains are already found in present CDS
#                min_number = int(rule.partition(",")[0])
#                required_matches = rule.partition("[")[2].partition("]")[0].split(",")
#                essential_matches = rule.rpartition("[")[2].partition("]")[0].split(",")
#                if essential_matches == ['']:
#                    essential_matches = []
#                cluster_results = list(set(required_matches) & set(cdsresults))
#                missing_results = min_number - len(cluster_results)
#                #If more than one, search nearby CDSs for the other needed domains
#                nrcds = 0
#                if len(cluster_results) > 0:
#                    locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
#                    nrcds = 1
#                    for othercds in [key for key in results_by_id.keys() if key != cds]:
#                        needed_hits_found = list(set([res.query_id for res in results_by_id[othercds]]) & set(required_matches))
#                        if len(needed_hits_found) > 0:
#                            feature = feature_by_id[othercds]
#                            flocations = [feature.location.start, feature.location.end]
#                            #If hit found in nearby CDS, add it to the set of found domains relevant to the present rule
#                            if min([abs(max(locations) - min(flocations)), abs(max(locations) - min(flocations))]) < cutoff:
#                                nrcds += 1
#                                cluster_results.extend(needed_hits_found)
#                                missing_results = missing_results - 1
#                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
#                if nrcds >= min_number and len(set(essential_matches) & set(cluster_results)) == len(essential_matches):
#                    if not (_type != "none" and clustertype == "other"):
#                        if _type == "none" or _type == "other" or _type == clustertype:
#                            _type = clustertype
#                        elif clustertype not in _type:
#                            _type = clustertype + "-" + _type
#                        break
#        #Save type to typedict
#        typedict[cds] = _type
#    return typedict


#def get_nseq():
#    nseqdict = {}
#    # FIXME: Implement this again
#    return nseqdict
#    for hmm in _signature_profiles:
#        hmmfile = hmm.hmm_file
#        for line in open(hmmfile, 'r'):
#            if line.startswith('NSEQ '):
#                nseqdict[hmm.name] = line[6:].strip()
#                break
#        if not hmm.name in nseqdict:
#            nseqdict[hmm.name] = "?"
#
#    return nseqdict


#def overlaps(feature1, feature2):
#    if (feature2.location.start <= feature1.location.start <= feature2.location.end) or (feature2.location.start <= feature1.location.end <= feature2.location.end):
#        return True
#    else:
#        return False


#def remove_irrelevant_allorfs(seq_record):
#    #Get features
#    allfeatures = utils.get_cds_features(seq_record)
#    #Remove auto-orf features without unique sec_met qualifiers; remove glimmer ORFs overlapping with sec_met auto-orfs not catched by Glimmer
#    auto_orf_features = [feature for feature in allfeatures if 'note' in feature.qualifiers and "auto-all-orf" in feature.qualifiers['note']]
#    other_features = [feature for feature in allfeatures if 'note' not in feature.qualifiers or "auto-all-orf" not in feature.qualifiers['note']]
#    to_delete = []
#    for autofeature in auto_orf_features:
#        if 'sec_met' not in autofeature.qualifiers:
#            to_delete.append(autofeature)
#        else:
#            glimmer_has_sec_met = False
#            for otherfeature in other_features:
#                if overlaps(autofeature, otherfeature) and 'sec_met' in otherfeature.qualifiers:
#                    to_delete.append(autofeature)
#                    glimmer_has_sec_met = True
#            if glimmer_has_sec_met is False:
#                for otherfeature in other_features:
#                    if overlaps(autofeature, otherfeature) and 'secment' not in otherfeature.qualifiers:
#                        to_delete.append(otherfeature)
#    featurenrs = []
#    idx = 0
#    for feature in seq_record.features:
#        if feature in to_delete:
#            featurenrs.append(idx)
#        idx += 1
#    featurenrs.reverse()
#    for featurenr in featurenrs:
#        del seq_record.features[featurenr]


#def add_additional_nrpspks_genes(typedict, results_by_id, seq_record, nseqdict):
#    '''Add annotations for additional NRPS/PKS genes in the cluster'''
#    nrpspksdomains = [
#        "PKS_KS",
#        "PKS_AT",
#        "ATd",
#        "ene_KS",
#        "mod_KS",
#        "hyb_KS",
#        "itr_KS",
#        "tra_KS",
#        "Condensation",
#        "AMP-binding",
#        "A-OX"
#    ]
#    clustercdsfeatures = utils.get_withincluster_cds_features(seq_record)
#    othercds_with_results = [
#        cds for cds in clustercdsfeatures
#        if utils.get_gene_id(cds) in results_by_id and typedict[utils.get_gene_id(cds)] == "none"
#    ]
#    for cds in othercds_with_results:
#        cdsresults = [res.query_id for res in results_by_id[utils.get_gene_id(cds)]]
#        if len(set(nrpspksdomains) & set(cdsresults)) >= 1:
#            _update_sec_met_entry(cds, results_by_id[utils.get_gene_id(cds)], "other", nseqdict)

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
