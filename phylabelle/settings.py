"""
.. _settings-detailed:

Settings
========

Note that all specified relative paths should relate to your projects root
directory.

.. _settings-phylophlan:

.. describe:: PHYLOPHLAN

    Set this to the path, where PhyloPhlAn is installed to.

    .. code-block:: python

        PHYLOPHLAN = '/home/USER/PATH/TO/nsegata-phylophlan-f2d78771d71d/'


.. describe:: STD_LINE_LEN

    length of a line for help text output

.. describe:: EMAIL

    ftp download needs an email address as identifier

.. describe:: DIRECTORIES

    this dictionary describes the folder-structure, used for a phylabelle
    project. It is possible to adapt every entry to your liking

.. describe:: DB

    phylabelle uses an sqlite database for internal storage. This
    variable states the name of the database that is used.

.. describe:: REFSEQ_TABLE, GENBANK_TABLE

    These values are the names of the local copies that phylabelle makes of the
    summary tables of refseq/genbank.

.. describe:: NAME_ILLEGALS

    In the name mapping process organism names are stripped from certain
    very common substrings, that are not unique and not necessary in an
    organism-name. Without these illegal strings, phylabelle would not notice,
    that, e.g. escherichia coli K12 and escherichia coli str. K12 refer to the
    same organism.

.. describe:: FTP_REPOSITORY

    base url of the FTP proteome source

.. describe:: FTP_SOURCES

    this dict enlists the sources of the summary files, that enable the
    determination of the exact location of a proteome.

.. describe:: LS_ASM_HEADER

    this describes, which fields of an Assembly should be displayed on request.

"""

import logging
import os
from collections import OrderedDict

#: location of PhyloPhlAn
PHYLOPHLAN = '/home/USER/PATH/TO/nsegata-phylophlan-f2d78771d71d/'

#: length of a line for help text output
STD_LINE_LEN = 80

EMAIL = 'user@host.whatever'

DIRECTORIES = {
    'proteomes': 'data/proteomes',
    'proteome_seq': 'data/proteomes/fasta',
    # directory to move the finished phylogenetic tree to
    'phylo': 'data/phylo',
}

#: phylabelle uses an sqlite filebased db
DB = 'data/data.db'

#: in the name mapping process organism names are stripped from certain
#: very common substrings, that are not unique and not necessary in an organism-name.
#: These are listed here.
NAME_ILLEGALS = ['sp.', 'subsp.', 'str.', 'serovar', 'strain=']

REFSEQ_TABLE = 'assembly_summary_refseq.txt'
GENBANK_TABLE = 'assembly_summary_genbank.txt'


#: main repository for proteome downloads
FTP_REPOSITORY = 'ftp.ncbi.nih.gov'

#: use paths relative to FTP_REPOSITORY
FTP_SOURCES = {
    REFSEQ_TABLE: '/genomes/refseq/assembly_summary_refseq.txt',
    GENBANK_TABLE: '/genomes/genbank/assembly_summary_genbank.txt'
}

#: Header for Assembly list output, mapped to
#: attributes of the Assembly object
LS_ASM_HEADER = OrderedDict([
    ('Accession', 'accession'),
    ('Name', 'organism_name'),
    ('Strain', 'infraspecific_name'),
    ('Label', 'label'),
    ('Assembly Level', 'assembly_level')
])

# logger settings
VERBOSE_LEVEL = logging.INFO
LOGFILE = 'main.log'
LOGFILE = os.path.join(os.getcwd(), LOGFILE)

#: indicates ranking between assembly levels
ASSEMBLY_LEVELS = ['Contig', 'Scaffold', 'Chromosome with gaps', 'Chromosome',
                   'Complete Genome']
