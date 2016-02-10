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

import helperlibs
import Bio

name = "local_file"
short_description = "{}: Load sequences from a local file".format(name)
priority = 1


def check_prereqs(options):
    '''Check for prerequisites'''
    failure_messages = []

    return failure_messages


def get_versions(options):
    '''Get all utility versions'''
    versions = [
        'Biopython {}'.format(Bio.__version__),
        'helperlibs {}'.format(helperlibs.version)
    ]
    return versions
