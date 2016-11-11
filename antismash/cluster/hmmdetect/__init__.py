# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2016 Kai Blin
# Danish Technical University
# Novo Nordisk Foundation Center for Biosustainability
# Novel Bioactive Compounds
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Load sequences from a local file

"""
import logging
import Bio
from Bio.SeqFeature import SeqFeature, FeatureLocation
from antismash import (
    config,
    utils,
)

name = "hmmdetect"
short_description = "{}: Detect secondary metabolite gene cluster".format(name)
priority = 1

_required_binaries = [
    ('hmmsearch', '3.1b1'),
    ('hmmpress', '3.1b1'),
]

_required_files = [
    'hmmdetails.txt',
    'cluster_rules.txt',
    'filterhmmdetails.txt',
]

_markov_model = [
    'bgc_seeds.hmm',
]

_binary_extensions = [
    '.h3f',
    '.h3i',
    '.h3m',
    '.h3p',
]


class HmmSignature(object):
    """HMM signature"""
    def __init__(self, name, description, cutoff, hmm_file):
        self.hmm_file = utils.get_full_path(__file__, hmm_file)
        self.name = name
        self.description = description
        self.cutoff = cutoff

#Define all profiles using details from hmmdetails.txt
hmmdetails = [
    _line.split("\t")
    for _line in open(utils.get_full_path(__file__, "hmmdetails.txt"), "r").read().split("\n")
    if _line.count("\t") == 3
]
_signature_profiles = [
    HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails
]


def check_prereqs(options):
    '''Check for prerequisites'''
    failure_messages = []
    for binary_name, binary_version in _required_binaries:
        if utils.locate_executable(binary_name) is None:
            failure_messages.append('Failed to locate executable for {!r}'.format(binary_name))
        # TODO: Check binary version here

    for hmm in _markov_model:
        hmm = utils.get_full_path(__file__, hmm)
        if utils.locate_file(hmm) is None:
            failure_messages.append('Failed to locate file {!r}'.format(hmm))
            continue
        for ext in _binary_extensions:
            binary = "{}{}".format(hmm, ext)
            if utils.locate_file(binary) is None:
                _, err, retcode = utils.run_hmmpress(hmm)
                if retcode != 0:
                    failure_messages.append('Failed to hmmpress {!r}: {!r}'.format(hmm, err))
                break

    for file_name in _required_files:
        full_name = utils.get_full_path(__file__, file_name)
        if utils.locate_file(full_name) is None:
            failure_messages.append('Failed to locate file {!r}'.format(full_name))


    lineno = 1
    for line in open(utils.get_full_path(__file__, 'cluster_rules.txt'), 'r'):
        if line.count('\t') != 3:
            failure_messages.append('Failed to use cluster rules from line {lineno} due to misformatting:\n{line!r}'.format(
                line=line, lineno=lineno))
        lineno += 1

    lineno = 1
    for line in open(utils.get_full_path(__file__, 'hmmdetails.txt'), 'r'):
        if line.count('\t') != 3:
            failure_messages.append('Failed to load HMM profile settings from line {lineno} due to misformatting:\n{line!r}'.format(
                line=line, lineno=lineno))
        lineno += 1


    return failure_messages


def detect(seq_record, options):
    '''Detect signature genes using HMM profiles'''
    logging.info('Detecting gene clusters using HMM library')
    feature_by_id = utils.get_feature_dict(seq_record)
    full_fasta = utils.get_multifasta(seq_record)
    rulesdict = create_rules_dict()
    results = []
    sig_by_name = {}
    results_by_id = {}
    runresults = utils.run_hmmsearch(utils.get_full_path(__file__, 'bgc_seeds.hmm'), full_fasta, use_tempfile=True)
    for sig in _signature_profiles:
        sig_by_name[sig.name] = sig

    logging.debug(sig_by_name.keys())
    for runresult in runresults:
        #Store result if it is above cut-off
        for hsp in runresult.hsps:
            try:
                sig = sig_by_name[hsp.query_id]
            except KeyError:
                logging.error('BUG: Failed to find signature for ID %s', hsp.query_id)
                continue
            if hsp.bitscore > sig.cutoff:
                results.append(hsp)
                if hsp.hit_id not in results_by_id:
                    results_by_id[hsp.hit_id] = [hsp]
                else:
                    results_by_id[hsp.hit_id].append(hsp)

    #Filter results by comparing scores of different models (for PKS systems)
    results, results_by_id = filter_results(results, results_by_id)

    #Use rules to determine gene clusters
    typedict = apply_cluster_rules(results_by_id, feature_by_id, rulesdict, inclusive=options.inclusive)

    #Find number of sequences on which each pHMM is based
    nseqdict = get_nseq()

    #Save final results to seq_record
    for cds in results_by_id.keys():
        feature = feature_by_id[cds]
        if typedict[cds] != "none":
            _update_sec_met_entry(feature, results_by_id[cds], typedict[cds], nseqdict)

    find_clusters(seq_record, rulesdict)

    #Find additional NRPS/PKS genes in gene clusters
    add_additional_nrpspks_genes(typedict, results_by_id, seq_record, nseqdict)

    #Add details of gene cluster detection to cluster features
    store_detection_details(rulesdict, seq_record)

    #If all-orfs option on, remove irrelevant short orfs
    if options.all_orfs:
        remove_irrelevant_allorfs(seq_record)


def get_versions(options):
    '''Get all utility versions'''
    # TODO actually calculate this
    hmmer_version = '3.1b1'
    hmmer_multithreading = True
    versions = [
        'Biopython {}'.format(Bio.__version__),
        'hmmsearch {version} (multithreading enabled: {multithreading})'.format(
            version=hmmer_version, multithreading=hmmer_multithreading),
    ]
    return versions


def find_clusters(seq_record, rulesdict):
    #Functions that detects the gene clusters based on the identified core genes
    features = utils.get_cds_features(seq_record)
    clustertype = ""
    clusters = []
    cfg = config.get_config()
    clusternr = cfg.next_clusternr
    state = "seed"
    last_cluster_extended = False

    for feature in features:
        if state == "seed":
            if 'sec_met' not in feature.qualifiers or len([feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat]) == 0:
                continue
            coregenetype = [feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat][0].partition("Type: ")[2]
            clustertype = coregenetype
            if "-" in clustertype:
                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
                extension = max([rulesdict[value][2] for value in clustertype.split("-")])
            else:
                cutoff = rulesdict[clustertype][1]
                extension = rulesdict[clustertype][2]
            start = max(feature.location.start - extension, 0)
            if int(start) < int(feature.location.end):
                loc = FeatureLocation(start, feature.location.end)
            else:
                loc = FeatureLocation(feature.location.end, start)
            new_cluster = SeqFeature(loc, type="cluster")
            new_cluster.qualifiers['product'] = [clustertype]
            new_cluster.qualifiers['note'] = ["Cluster number: " + str(clusternr)]
            clusternr += 1
            clusters.append(new_cluster)
            state = "extend"
            last_cluster_extended = False

        elif state == "extend":
            cluster = clusters[-1]
            if "-" in clustertype:
                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
                extension = max([rulesdict[value][2] for value in clustertype.split("-")])
            else:
                cutoff = rulesdict[clustertype][1]
                extension = rulesdict[clustertype][2]
            if 'sec_met' not in feature.qualifiers:
                if feature.location.start > cluster.location.end + cutoff:
                    # ok, no hits for too long, done with this cluster
                    state = "seed"
                    end = min(len(seq_record), cluster.location.end + extension)
                    if int(cluster.location.start) < int(end):
                        cluster.location = FeatureLocation(cluster.location.start, end)
                    else:
                        cluster.location = FeatureLocation(end, cluster.location.start)
                    last_cluster_extended = True
                continue
            coregenetype = [feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat][0].partition("Type: ")[2]
            if coregenetype not in clustertype and coregenetype != 'other':
                if clustertype == "other" and coregenetype != "other":
                    clustertype = coregenetype
                else:
                    if "-" not in coregenetype:
                        if coregenetype not in clustertype and coregenetype != "other":
                            clustertype += "-" + coregenetype
                    else:
                        for partctype in [ctype for ctype in coregenetype.split("-") if ctype != "other"]:
                            if partctype not in clustertype:
                                clustertype += "-" + partctype
                cluster.qualifiers['product'] = [clustertype]
            if int(cluster.location.start) < int(feature.location.end):
                cluster.location = FeatureLocation(cluster.location.start, feature.location.end)
            else:
                cluster.location = FeatureLocation(feature.location.end, cluster.location.start)

    if len(clusters) > 0 and not last_cluster_extended:
        cluster = clusters[-1]
        end = min(len(seq_record), feature.location.end + extension)
        if int(cluster.location.start) < int(end):
            cluster.location = FeatureLocation(cluster.location.start, end)
        else:
            cluster.location = FeatureLocation(end, cluster.location.start)

    seq_record.features.extend(clusters)
    cfg.next_clusternr = clusternr


def filter_results(results, results_by_id):
    #Filter results by comparing scores of different models (for PKS systems)
    for line in open(utils.get_full_path(__file__, "filterhmmdetails.txt"), "r").read().split("\n"):
        filterhmms = line.split(",")
        for cds in results_by_id.keys():
            cdsresults = results_by_id[cds]
            hmmhits = [hit.query_id for hit in cdsresults]
            #Check if multiple competing HMM hits are present
            competing_hits = set(hmmhits) & set(filterhmms)
            if len(competing_hits) > 1:
                #Identify overlapping hits
                overlapping_groups = []
                for hit in cdsresults:
                    for otherhit in [cdsresult for cdsresult in cdsresults if hit != cdsresult]:
                        overlap = len(set(range(hit.hit_start, hit.hit_end)) & set(range(otherhit.hit_start, otherhit.hit_end)))
                        if overlap > 20:
                            added = "n"
                            for group in overlapping_groups:
                                if hit in group and otherhit in group:
                                    added = "y"
                                    break
                                elif hit in group and otherhit not in group:
                                    group.append(otherhit)
                                    added = "y"
                                    break
                                elif hit not in group and otherhit in group:
                                    group.append(hit)
                                    added = "y"
                                    break
                            if added == "n":
                                overlapping_groups.append([hit, otherhit])
                #Remove worst-scoring of overlapping hits
                for group in overlapping_groups:
                    highestscore = max([hit.bitscore for hit in group])
                    hit_with_highestscore = group[[hit.bitscore for hit in group].index(highestscore)]
                    to_delete = [hit for hit in group if hit != hit_with_highestscore]
                    for res in [res for res in results]:
                        if res in to_delete:
                            del results[results.index(res)]
                            del results_by_id[cds][results_by_id[cds].index(res)]
    return results, results_by_id


def create_rules_dict():
    "Create a cluster rules dictionary from the cluster rules file"
    rulesdict = {}
    first = True
    #TODO: We should move all user-customizable files into config subdirectory; the rulefiles are redundant also in hmm_detection_dblookup
    for line in open(utils.get_full_path(__file__, "cluster_rules.txt"), "r"):
        # skip the first line with the legend
        if first:
            first = False
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        key = parts.pop(0)
        rules = parts.pop(0)
        cutoff = int(parts.pop(0)) * 1000
        extension = int(parts.pop(0)) * 1000
        rulesdict[key] = (rules, cutoff, extension)
    return rulesdict


def apply_cluster_rules(results_by_id, feature_by_id, rulesdict, inclusive):
    "Apply cluster rules to determine if HMMs lead to secondary metabolite core gene detection"
    typedict = {}
    for cds in results_by_id.keys():
        _type = "none"
        cdsresults = [res.query_id for res in results_by_id[cds]]
        for clustertype in rulesdict.keys():
            if not inclusive and clustertype.startswith('cf_'):
                continue
            single_rules = [rule for rule in rulesdict[clustertype][0].split(" or ") if " & " not in rule and "cluster(" not in rule]
            combined_rules = [rule[1:-1] for rule in rulesdict[clustertype][0].split(" or ") if " & " in rule]
            cluster_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "cluster(" in rule]
            minimum_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "minimum(" in rule]
            if "-" in clustertype:
                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
            else:
                cutoff = rulesdict[clustertype][1]
            #Assign cluster type if a single argument rule matches
            #Example rule format: "Domain1"
            if len(set(single_rules) & set(cdsresults)) >= 1:
                if not (_type != "none" and clustertype == "other"):
                    if _type == "none" or _type == "other" or _type == clustertype:
                        _type = clustertype
                    elif clustertype not in _type:
                        _type = clustertype + "-" + _type
                if _type != "other":
                    continue
            #Assign cluster type if a combinatorial argument rule matches
            #Example rule format: "(Domain1 & Domain2)"
            for rule in combined_rules:
                required_matches = rule.split(" & ")
                if len(set(required_matches) & set(cdsresults)) == len(required_matches):
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type:
                            _type = clustertype + "-" + _type
            if _type == clustertype and _type != "other":
                continue
            #Assign cluster type if distance-based combinatorial parameter matches
            #Example rule format: "cluster(Domain1,Domain2)"
            for rule in cluster_rules:
                #Find which needed domains are already found in present CDS
                required_matches = rule.split(",")
                cluster_results = list(set(required_matches) & set(cdsresults))
                missing_results = list(set(required_matches) - set(cdsresults))
                #If more than one, search nearby CDSs for the other needed domains
                if len(cluster_results) > 0:
                    locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
                    for othercds in results_by_id.keys():
                        needed_hits_found = list(set([res.query_id for res in results_by_id[othercds]]) & set(missing_results))
                        if len(needed_hits_found) > 0:
                            feature = feature_by_id[othercds]
                            flocations = [feature.location.start, feature.location.end]
                            #If hit found in nearby CDS, add it to the set of found domains relevant to the present rule
                            if min([abs(max(locations) - min(flocations)), abs(max(locations) - min(flocations))]) < cutoff:
                                cluster_results.extend(needed_hits_found)
                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
                if len(list(set(required_matches) - set(cluster_results))) == 0:
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type:
                            _type = clustertype + "-" + _type
                        break
            #Assign cluster type if distance-based combinatorial parameter matches with minimum number of domain hits
            #Example rule format: "minimum(4,[Domain1,Domain2], [Domain1])"
            #This denotes that at least 4 hits of either Domain1 or Domain2 should be found in the same region with at least 1 hit from Domain 1
            for rule in minimum_rules:
                #Find which needed domains are already found in present CDS
                min_number = int(rule.partition(",")[0])
                required_matches = rule.partition("[")[2].partition("]")[0].split(",")
                essential_matches = rule.rpartition("[")[2].partition("]")[0].split(",")
                if essential_matches == ['']:
                    essential_matches = []
                cluster_results = list(set(required_matches) & set(cdsresults))
                missing_results = min_number - len(cluster_results)
                #If more than one, search nearby CDSs for the other needed domains
                nrcds = 0
                if len(cluster_results) > 0:
                    locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
                    nrcds = 1
                    for othercds in [key for key in results_by_id.keys() if key != cds]:
                        needed_hits_found = list(set([res.query_id for res in results_by_id[othercds]]) & set(required_matches))
                        if len(needed_hits_found) > 0:
                            feature = feature_by_id[othercds]
                            flocations = [feature.location.start, feature.location.end]
                            #If hit found in nearby CDS, add it to the set of found domains relevant to the present rule
                            if min([abs(max(locations) - min(flocations)), abs(max(locations) - min(flocations))]) < cutoff:
                                nrcds += 1
                                cluster_results.extend(needed_hits_found)
                                missing_results = missing_results - 1
                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
                if nrcds >= min_number and len(set(essential_matches) & set(cluster_results)) == len(essential_matches):
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type:
                            _type = clustertype + "-" + _type
                        break
        #Save type to typedict
        typedict[cds] = _type
    return typedict


def get_nseq():
    nseqdict = {}
    # FIXME: Implement this again
    return nseqdict
    for hmm in _signature_profiles:
        hmmfile = hmm.hmm_file
        for line in open(hmmfile, 'r'):
            if line.startswith('NSEQ '):
                nseqdict[hmm.name] = line[6:].strip()
                break
        if hmm.name not in nseqdict:
            nseqdict[hmm.name] = "?"

    return nseqdict


def overlaps(feature1, feature2):
    if (feature2.location.start <= feature1.location.start <= feature2.location.end
            or feature2.location.start <= feature1.location.end <= feature2.location.end):
        return True
    else:
        return False


def remove_irrelevant_allorfs(seq_record):
    #Get features
    allfeatures = utils.get_cds_features(seq_record)
    #Remove auto-orf features without unique sec_met qualifiers; remove glimmer ORFs overlapping with sec_met auto-orfs not catched by Glimmer
    auto_orf_features = [feature for feature in allfeatures if 'note' in feature.qualifiers and "auto-all-orf" in feature.qualifiers['note']]
    other_features = [feature for feature in allfeatures if 'note' not in feature.qualifiers or "auto-all-orf" not in feature.qualifiers['note']]
    to_delete = []
    for autofeature in auto_orf_features:
        if 'sec_met' not in autofeature.qualifiers:
            to_delete.append(autofeature)
        else:
            glimmer_has_sec_met = False
            for otherfeature in other_features:
                if overlaps(autofeature, otherfeature) and 'sec_met' in otherfeature.qualifiers:
                    to_delete.append(autofeature)
                    glimmer_has_sec_met = True
            if glimmer_has_sec_met is False:
                for otherfeature in other_features:
                    if overlaps(autofeature, otherfeature) and 'secment' not in otherfeature.qualifiers:
                        to_delete.append(otherfeature)
    featurenrs = []
    idx = 0
    for feature in seq_record.features:
        if feature in to_delete:
            featurenrs.append(idx)
        idx += 1
    featurenrs.reverse()
    for featurenr in featurenrs:
        del seq_record.features[featurenr]


def add_additional_nrpspks_genes(typedict, results_by_id, seq_record, nseqdict):
    '''Add annotations for additional NRPS/PKS genes in the cluster'''
    nrpspksdomains = [
        "PKS_KS",
        "PKS_AT",
        "ATd",
        "ene_KS",
        "mod_KS",
        "hyb_KS",
        "itr_KS",
        "tra_KS",
        "Condensation",
        "AMP-binding",
        "A-OX"
    ]
    clustercdsfeatures = utils.get_withincluster_cds_features(seq_record)
    othercds_with_results = [
        cds for cds in clustercdsfeatures
        if utils.get_gene_id(cds) in results_by_id and typedict[utils.get_gene_id(cds)] == "none"
    ]
    for cds in othercds_with_results:
        cdsresults = [res.query_id for res in results_by_id[utils.get_gene_id(cds)]]
        if len(set(nrpspksdomains) & set(cdsresults)) >= 1:
            _update_sec_met_entry(cds, results_by_id[utils.get_gene_id(cds)], "other", nseqdict)


def store_detection_details(rulesdict, seq_record):
    '''Store the details about why a cluster was detected'''
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        type_combo = utils.get_cluster_type(cluster)
        if '-' in type_combo:
            clustertypes = type_combo.split('-')
        else:
            clustertypes = [type_combo]

        if 'note' not in cluster.qualifiers:
            cluster.qualifiers['note'] = []
        rule_string = "Detection rule(s) for this cluster type:"
        for clustertype in clustertypes:
            rule_string += " %s: (%s);" % (clustertype, rulesdict[clustertype][0])

        cluster.qualifiers['note'].append(rule_string)


def _update_sec_met_entry(feature, results, clustertype, nseqdict):
    '''Update the sec_met entry for a feature'''
    result = "; ".join(["%s (E-value: %s, bitscore: %s, seeds: %s)" % (
        res.query_id, res.evalue, res.bitscore, nseqdict.get(res.query_id, '?')) for res in results])

    if 'sec_met' not in feature.qualifiers:
        feature.qualifiers['sec_met'] = [
            "Type: %s" % clustertype,
            "Domains detected: %s" % (result),
            "Kind: biosynthetic"
        ]
    else:
        for ann in feature.qualifiers['sec_met']:
            if not ann.startswith("Domains detected"):
                continue
            ann += "Domains detected: %s" % (result)
