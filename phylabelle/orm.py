"""
ORM
===

Database stuff.
"""
from sqlalchemy import Table, Column, String, ForeignKey, Date, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import mapper, relationship

from phylabelle.fileio import format_get_lines

# store mappings in dict:
mappings = {}

Base = declarative_base()

assemblies = Table('assemblies', Base.metadata,
                   # bioproject id
                   Column('accession', String(20), primary_key=True),
                   Column('label', Boolean),
                   Column('species_name', String(40)),
                   Column('bioproject_labelsource', String(20)),
                   Column('bioproject', String(20)),
                   Column('assembly_name', String(20)),
                   Column('organism_name', String(100)),
                   Column('infraspecific_name', String(100)),
                   Column('version_status', String(20)),
                   # this is the Taxon-ID of the strain, usually not a foreign key
                   Column('tax_id', String(20)),
                   Column('species_tax_id', String(20)),
                   Column('assembly_level', String(20)),
                   Column('ftp_source', String(100)),
                   Column('assembly_level', String(30)),
                   Column('submission_date', Date),
                   Column('seq_file', String(50)),
                   Column('report', String(50)),
                   Column('feature_table', String(50))
                   )

#: store those assemblies which have been left out, as reference to
#: chosen assembly
rejected_assemblies = Table('other_assemblies', Base.metadata,
                            Column('accession', String(20), primary_key=True),
                            Column('ftp_source', String(20)),
                            Column('_assembly_acc', String(20), ForeignKey('assemblies.accession'))
                            )

nucleotides = Table('nucleotides', Base.metadata,
                    # ID given
                    Column('accession', String(20), primary_key=True),
                    Column('type', String(20)),
                    Column('_assembly_acc', String(20), ForeignKey('assemblies.accession'))
                    )


class Assembly(object):
    """
    Basic organism entity
    """

    #: provide a short description which values the attributes
    # are supposed to contain
    _help = {
        'accession': 'assembly accession',
        'bioproject': 'bioproject id of primary submission. In case of Refseq assemblies, \
                       the field still contains the bioproject id of the primary submission \
                       linked, since these were used to assign the labels',
        'assembly_name': 'assembly name',
        'tax_id': 'Most specific NCBI Taxonomy ID',
        'species_tax_id': 'NCBI Taxonomy ID on species level',
        'organism_name': 'organism name. Note that the values in this field do not necessarily \
                          contain a strain identifier, but some do.',
        'infraspecific_name': 'infraspecific name, like strain identifier. Does not \
                               necessarily contain a value',
    }

    @classmethod
    def from_proto(cls, proto_asm, metadata, files):
        """
        generate Assembly object from proto_assembly and metadata
        :param ProtoAssembly proto_asm:
        :param dict metadata: metadata from assembly_summary_tables
        :param iterable files: iterable of filenames associated to assembly
        :return Assembly asm:
        """
        asm = Assembly(metadata['accession'])
        try:
            asm.label_source_bioproject = proto_asm.bioproj_acc
        except AttributeError:
            pass

        asm.label = proto_asm.label
        try:
            asm.species_name = proto_asm.species_name
        except AttributeError:
            pass

        asm.bioproject = metadata['bioproject']
        asm.assembly_name = metadata['assembly_name']
        asm.submission_data = metadata['submission_date']
        asm.version_status = metadata['version_status']
        asm.assembly_level = metadata['assembly_level']
        asm.ftp_source = metadata['ftp_source']
        asm.tax_id = metadata['tax_id']
        asm.species_tax_id = metadata['species_tax_id']
        asm.organism_name = metadata['organism_name']
        asm.infraspecific_name = metadata['infraspecific_name']
        asm.nucleotides = []
        asm.other = proto_asm.assemblies

        for file_name in files:
            if 'protein.faa' in file_name:
                asm.seq_file = file_name
            elif 'report' in file_name:
                asm.report = file_name
            else:
                asm.feature_table = file_name

        return asm

    @classmethod
    def iter_all(cls, session, subset=None):
        """
        iterate over all available assemblies
        :param session: db-session
        :param list subset: list of accessions; load only ids, contained in this list
        """
        if filter is not None:
            asms = session.query(cls).filter(cls.accession.in_(subset)).all()
        else:
            asms = session.query(cls).all()

        return iter(asms)

    @classmethod
    def get_species_index(cls, session):
        """
        get a dict containing accession as key and species_tax_id as value
        :return:
        """
        return dict(session.query(cls.accession, cls.species_tax_id).all())

    @property
    def refseq(self):
        """
        :return: boolean, indicating whether asm is a refseq assembly
        """
        return 'GCF' in self.accession

    @classmethod
    def show_fields(cls):
        """
        show a description of searchable fields
        """
        for field_name, help_text in cls._help.iteritems():
            it = format_get_lines(help_text)
            print '{:20}:\t{:10}'.format(field_name, it.next())

            for line in it:
                print '{:20}\t{:10}'.format('', line)

    def __init__(self, accession):
        """
        :param accession: assembly accession
        """
        self.accession = accession

    def get_nucleotides(self, type_):
        """
        :param type_: type of nucleotide
        """
        for nucl in self.nucleotides:
            if nucl['type'] == type_:
                yield nucl

    def count_nucleotides(self, type_):
        """
        :param str type_: the type of nucleotide which shall be counted
        """
        return len(list(self.get_nucleotides(type_)))

    @property
    def num_plasmids(self):
        """
        :return: number of plasmids
        """
        return self.count_nucleotides('Plasmid')


class Nucleotide(object):
    def __init__(self, accession, type):
        self.accession = accession
        self.type = type


class RejectedAssembly(object):
    """
    Holds information about assemblies which have not been chosen,
    in favor of better ones
    """
    def __init__(self, accession, ftp_source):
        self.accession = accession
        self.ftp_source = ftp_source


"""
Mappings
========
bind objects to tables, set relationships
"""

mappings['Assembly'] = mapper(Assembly, assemblies)
mappings['Nucleotide'] = mapper(Nucleotide, nucleotides)
mappings['RejectedAssembly'] = mapper(RejectedAssembly, rejected_assemblies)

# Now add relationships:
mappings['Assembly'].add_properties({
    'nucleotides': relationship(Nucleotide, backref='assembly'),
    'rejected_assemblies': relationship(RejectedAssembly, backref='assembly'),
})
