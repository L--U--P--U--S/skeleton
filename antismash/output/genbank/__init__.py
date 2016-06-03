'''GenBank output module'''
import logging
import warnings
from os import path

from helperlibs.bio import seqio
from antismash import utils

name = "genbank"
short_description = "GenBank output"
priority = 1


def write(seq_records, options):
    '''Write GenBank records'''
    basename = seq_records[0].id
    output_name = path.join(options.outputfoldername, "%s.final.gbk" % basename)
    logging.debug("Writing %s seq_records to %r", len(seq_records), output_name)
    seqio.write(seq_records, output_name, 'genbank')
    i = 1
    for rec in seq_records:
        # For compatibility with the database importer, we have to check whether we are dealing
        # with a seq_record obtained from a file (then its class will be SeqRecord) or from a#
        # database (then its class will be DBSeqRecord)
        #
        # running the cluster extraction on a DBSeqRecord will throw an exception, as splitting
        # the object is not supported
        if rec.__class__.__name__ == 'DBSeqRecord':
            continue
        for cluster in utils.get_cluster_features(rec):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cluster_rec = rec[cluster.location.start:cluster.location.end]
            output_name = path.join(options.outputfoldername,
                                    "%s.cluster%03d.gbk" % (basename, i))
            seqio.write([cluster_rec], output_name, 'genbank')
            i += 1
