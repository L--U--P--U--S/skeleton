'''SeqRecord access helpers'''
import logging
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from antismash.utils.errors import InvalidLocationError


def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for feature in seq_record.features:
        if feature.type in types:
            features.append(feature)
    return features


def get_all_features_of_type_with_query(seq_record, feature_type, query_tag, query_value):
    """Return all features of type 'type' which contain a 'query_tag' with value 'query_value'
    Note: query has to be exact!"""

    features_with_type = get_all_features_of_type(seq_record, feature_type)
    features = []
    for feature_to_test in features_with_type:

        if (query_tag in feature_to_test.qualifiers) and (query_value in feature_to_test.qualifiers[query_tag]):
            features.append(feature_to_test)
    return features


def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")


def get_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    return get_all_features_of_type(seq_record, "cluster")


def get_withincluster_cds_features(seq_record):
    '''Get all CDS features within a cluster feature'''
    features = get_cds_features(seq_record)
    clusters = get_cluster_features(seq_record)
    withinclusterfeatures = []
    for feature in features:
        for cluster in clusters:
            if not (cluster.location.start <= feature.location.start <= cluster.location.end or
                    cluster.location.start <= feature.location.end <= cluster.location.end):
                continue
            if feature not in withinclusterfeatures:
                withinclusterfeatures.append(feature)
    return withinclusterfeatures


def get_cluster_type(cluster):
    "Get product type of a gene cluster"
    return cluster.qualifiers['product'][0]


def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id


def get_multifasta(seq_record):
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = get_cds_features(seq_record)
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers['translation'][0]
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-", ""), generic_protein)

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping", gene_id)
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta


def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"


def get_gene_annotation(feature):
    "Get the gene annotation from the product qualifier"
    if 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    return "unannotated orf"


def get_gene_accession(feature):
    "Get the gene accession number from protein_id, locus_tag or gene name, in that order"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_tag_found"


def cmp_feature_location(first, second):
    "Compare two features by their start/end locations"
    ret = cmp(first.location.start, second.location.start)
    if ret != 0:
        return ret
    return cmp(first.location.end, second.location.end)


def sort_features(seq_record):
    "Sort features in a seq_record by their position"
    #Check if all features have a proper location assigned
    for feature in seq_record.features:
        if feature.location is None:
            if feature.id != "<unknown id>":
                logging.error("Feature '%s' has no proper location assigned", feature.id)
            elif "locus_tag" in feature.qualifiers:
                logging.error("Feature '%s' has no proper location assigned", feature.qualifiers["locus_tag"][0])
            else:
                logging.error("File contains feature without proper location assignment")
            raise InvalidLocationError(feature)
    #Sort features by location
    seq_record.features.sort(cmp=cmp_feature_location)


# shamelessly stolen from the "real" antismash
def features_overlap(a, b):
    "Check if two features have overlapping locations"
    astart = min([a.location.start, a.location.end])
    aend = max([a.location.start, a.location.end])
    bstart = min([b.location.start, b.location.end])
    bend = max([b.location.start, b.location.end])
    return (astart >= bstart and astart <= bend) or \
           (aend >= bstart and aend <= bend) or \
           (bstart >= astart and bstart <= aend) or \
           (bend >= astart and bend <= aend)
