'''Core antiSMASH functionality'''
import logging

from antismash.utils.access import (
    get_cluster_features,
    sort_features,
    get_all_features_of_type,
)


def run_antismash(seq_records, plugins, options):
    '''Run the antiSMASH pipeline'''
    options.record_idx = 1
    options.orig_record_idx = 1
    options.next_clusternr = 1
    #FIXME: properly implement this
    options.all_orfs = False
    options.inclusive = False
    for seq_record in seq_records:
        logging.info('Analysing record %d', options.record_idx)

        detect_gene_clusters(seq_record, plugins, options)

        # Cluster-specific and cluster-indepented analyses here
        sort_features(seq_record)

    store_results(seq_records, plugins, options)


def detect_gene_clusters(seq_record, plugins, options):
    '''Detect secondary metabolite gene clusters'''
    for plugin in plugins['cluster']:
        plugin.detect(seq_record, options)

    if options.debug:
        for feature in get_cluster_features(seq_record):
            logging.debug(feature)


def store_results(seq_records, plugins, options):
    '''Store all the results'''
    logging.debug('Writing the output for %s sequence records', len(seq_records))
    for plugin in plugins['output']:
        plugin.write(seq_records, options)
