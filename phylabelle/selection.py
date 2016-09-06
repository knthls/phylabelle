# coding: utf-8

import datetime
import os
import warnings
from collections import defaultdict

from sortedcontainers import SortedSet

from phylabelle.fileio import Table, Header, Line
from phylabelle.maintenance import update_indices
from phylabelle.utils import parse_bool, NoBooleanValueException


def get_indices(settings, update=False):
    """
    get absolute paths of current refseq and genbank assembly-lists
    :param module settings:
    :param bool update: if true, indices are updated from FTP repository
    :return: a tuple of paths to both files
    """
    if not update:
        try:
            rs_tab = os.path.join('data', settings.REFSEQ_TABLE)
            gb_tab = os.path.join('data', settings.GENBANK_TABLE)
            return gb_tab, rs_tab
        except IndexError:
            update = True

    if update:
        gb_tab, rs_tab = update_indices(settings)
        return gb_tab, rs_tab


def shrink_name(name, illegals):
    """
    Eliminates strings contained in illegals. Name is split into substrings,
    split positions are the whitespace positions, If a substring
    appears more than once, only one repeat is kept.
    :param str name:
    :param list illegals: list of illegal strings
    """
    if illegals is not None:
        name_parts = [x for x in name.split() if not x.lower() in illegals]
    else:
        name_parts = name.split()

    for i in xrange(len(name_parts) - 2):
        if name_parts[i].lower() == name_parts[i+1].lower():
            del name_parts[i]

    return ' '.join(name_parts)


def get_species_name(name):
    """
    extracts the species from a complete organism name
    :param str name: organism name
    :returns tuple: species name, suffix
    """
    pos = name.find(' ') + 1
    pos2 = name[pos:].find(' ') + pos
    if pos2 < pos:
        return name, None
    elif name[pos:pos2] == 'sp.':
        pos2 += 1
        pos3 = pos2 + name[pos2:].find(' ')
        if pos3 < pos2:
            return name, None

    return name[:pos2], name[pos2:]


def unify_name(species_name, strain, illegals=None):
    """
    remove ambiguous naming patterns, making exact name matching possible
    :param species_name:
    :param strain:
    :param illegals:
    :return:
    """
    spec_name, suffix = get_species_name(species_name)
    if suffix:
        suffix = shrink_name(suffix, illegals=illegals)
        if strain in suffix:
            strain = suffix

    if 'strain=' in strain:
        strain = strain[7:].strip()

    return spec_name, strain


def asm_score(submission_date, version_status, asm_level, refseq, assembly_levels):
    """
    calculate a score for an assembly, giving higher values to
        1. newer assemblies
        2. assemblies with better version status ()
    :param datetime.date submission_date:
    :param str version_status:
    :param str asm_level:
    :param bool refseq: indicates whether asm is from refseq or not
    :param list assembly_levels: list of assembly levels. Has to contain at least asm_level
    """
    assembly_level = assembly_levels.index(asm_level)
    try:
        val = submission_date.toordinal()
    except AttributeError:
        val = 0
    val += int(version_status == 'latest') * 10e6
    val += assembly_level * 10e5
    val += int(refseq) * 10e7
    return val


class InvalidFileFormatError(Exception):
    pass


class InvalidAssemblyException(Exception):
    pass


class SimpleMappingIndex(object):
    """
    Map bioprojects to their assemblies in order to find best one
    """
    def __init__(self):
        self.proto_assemblies = defaultdict(lambda: None)
        self.assemblies = []
        self.rs_refs = defaultdict(lambda: set([]))

    def add_proto_asm(self, proto_asm):
        """
        add a ProtoAssembly to both dicts
        :param TrivialProtoAssembly proto_asm:
        """
        self.proto_assemblies[proto_asm.accession] = proto_asm
        self.assemblies.append(proto_asm)

    def map(self, acc, score, rs_ref, metadata):
        """
        map assembly to ProtoAssembly instance
        :param acc: accession
        :param score: assembly score
        :param rs_ref: refseq accession
        :param metadata: metadata
        """
        p_asm = self.proto_assemblies[acc]
        if p_asm is not None:
            p_asm.add_asm(acc, score, metadata)


class MappingIndex(object):
    """
    Map bioprojects and organism names to assemblies, and choose best scoring assembly
    for each organism.
    Background is, that Bioprojects in general are no unique key for assemblies and
    can even contain multiple organisms.
    """
    def __init__(self):
        # holds a mapping to ProtoAssembly-Objects with bioproject-ids as key
        self.bp_key = defaultdict(lambda: None)
        # two layered mapping
        # key 1: species name
        # key 2: strain name
        self.name_key = defaultdict(lambda: defaultdict(lambda: []))
        # asm_mapping contains assembly_accessions and maps them to proto_assemblies
        self.asm_mapping = {}
        self.rs_refs = defaultdict(lambda: set([]))
        self.assemblies = []

    def add_proto_asm(self, proto_asm):
        """
        add a ProtoAssembly to both dicts
        :param ProtoAssembly proto_asm:
        """
        self.bp_key[proto_asm.bioproj_acc] = proto_asm
        self.name_key[proto_asm.species_name][proto_asm.strain].append(proto_asm)

        self.assemblies.append(proto_asm)

    def register_ref(self, key, asm):
        """
        register a refseq assembly accession linked to proto assembly
        :param str key: accession/reference to refseq
        :param ProtoAssembly asm:
        """
        self.rs_refs[key].add(asm)

    def iter_names(self):
        """
        provide an iterator name_key
        :return: generator over names, and values
        """
        for name, dict_ in self.name_key.iteritems():
            for inf_name, values in dict_.iteritems():
                yield name, inf_name, values

    def map(self, bioproject_id, species_name, strain, acc, score, rs_ref, metadata):
        """
        Map an assembly (from assembly index) to ProtoAssembly instance
        (which is supposed to represent exactly one organism)
        :param bioproject_id:
        :param species_name:
        :param strain:
        :param acc: assembly accession
        :param score: assembly score
        :param rs_ref: refseq assembly
        :param metadata:
        """
        proto_asm1 = self.bp_key[bioproject_id]
        proto_asm2 = self.name_key[species_name][strain]

        mapped = set([])

        if proto_asm1:
            mapped.add(proto_asm1)

        if len(proto_asm2) > 0:
            for asm in proto_asm2:
                if asm.has_strain():
                    mapped.add(asm)

        for asm in mapped:
            asm.add_asm(acc, score, metadata)
            self.register_ref(rs_ref, asm)

        mapped = list(mapped)
        return mapped

    def get_mapped(self):
        """
        iterate over all proto assemblies which hold a mapped assembly
        """
        for asm in self.assemblies:
            if len(asm.assemblies) > 0:
                yield asm

    def check_integrity(self):
        """
        Check whether mapping process has produced contradicting labels
        """
        with Table('contradictions.tsv', 'w',
                   header=Header(sep='\t',
                                 list_=['Name',
                                        'Strain',
                                        'AsmSpeciesNamePlain',
                                        'AsmStrainPlain',
                                        'AsmSpeciesNameUnif',
                                        'AsmStrainPlain',
                                        'Asm.Label'
                                        ])) as tab:
            iter_ = self.iter_names()

            for name, inf_name, asms in iter_:
                label = None
                for asm in asms:
                    if label is None:
                        label = asm.label
                    else:
                        if label != asm.label:
                            for asm_ in asms:
                                tab.write([name, inf_name, asm_.plain_names[0],
                                           asm_.plain_names[1], asm_.species_name,
                                           asm_.strain, asm_.label])
                            break


class KeyTuple(tuple):
    """
    Items in ProtoAssembly.assemblies must be hashable, so this is a class,
    exteding tuples to hash only their 2nd value (accession)
    """
    def __init__(self, *args, **kwargs):
        super(KeyTuple, self).__init__(*args, **kwargs)

    def __hash__(self):
        return hash(self[1])


class TrivialProtoAssembly(object):
    """
    Wraps a SortedSet of assemblies.
    One ProtoAssembly is supposed to represent one distinct organism
    """
    def __init__(self, acc, label):
        self.accession = acc
        assert type(label) == bool
        self.label = label

        #: supposed to contain tuples of the form
        #: (SCORE, ASSEMBLY_ACCESSION, METADATA)
        #: sorted from highest to lowest score
        self.assemblies = SortedSet(key=lambda x: -x[0])

    def add_asm(self, acc, score, metadata):
        """
        add reference to actual genbank/refseq assembly, including ftp-source
        :param str acc: accession
        :param int score: assembly-score, indicating which one is best
        :param dict metadata: optional metadata
        :return bool: indicating whether asm has been accepted (i.e.) not already
        contained in the references-set
        """
        old_size = len(self.assemblies)
        self.assemblies.add(KeyTuple((score, acc, metadata)))
        accepted = old_size < len(self.assemblies)
        return accepted


class ProtoAssembly(TrivialProtoAssembly):
    """
    Represents one distinct organism and holds references to all assemblies linked to it.
    During the download Process it yields the best scoring assemblies until a valid one
    (i.e. one with an associated, complete proteome) occurs.
    """
    def __init__(self, bioproj_acc, species_name, strain, label, name_illegals=None):
        super(ProtoAssembly, self).__init__(None, label)

        self.bioproj_acc = bioproj_acc
        if not len(strain) > 0:
            # print bioproj_acc, species_name, strain
            pass

        self.species_name, self.strain = unify_name(species_name, strain, illegals=name_illegals)
        self.plain_names = (species_name, strain)

    def has_strain(self):
        return len(self.strain) > 1

    def __eq__(self, other):
        try:
            return self.match_name(other.species_name, other.strain)
        except:
            return False

    def match_name(self, species_name, strain):
        """
        given a species name and a strain, determine whether these represent the same
        organism
        :param species_name:
        :param strain:
        :return:
        """
        m1 = self.species_name.lower().strip() == species_name.lower().strip()
        m2 = self.strain.lower().strip() in strain.lower().strip()
        m3 = strain.lower().strip() in self.strain.lower().strip()
        return m1 and (m2 or m3)


class PreProcessor(object):
    """
    Manager class for mapping names/labels/organisms to assemblies.
    """
    def __init__(self, label_file, mapping_log, settings, update=False):
        """
        :param str label_file: absolute path to file, containing 2 tabs Assembly-Accession and labels
        :param bool update: update assembly summaries
        """
        self.mapping_index = None
        self.read_labels(label_file)
        self.mapping_log = mapping_log
        self.gb_table, self.rs_table = get_indices(settings, update=update)
        self.assembly_levels = settings.ASSEMBLY_LEVELS
        self.name_illegals = settings.NAME_ILLEGALS

    def _read_line(self, line, valid_accessions=None):
        acc = line['assembly_accession']
        if valid_accessions is not None:
            if acc not in valid_accessions:
                raise InvalidAssemblyException

        plain_spec = line['organism_name']
        plain_strain = line['infraspecific_name']

        spec, strain = unify_name(plain_spec, plain_strain,
                                  illegals=self.name_illegals)

        bioproj = line['bioproject']

        rs_ref = line['gbrs_paired_asm']

        version_status = line['version_status'].strip()
        asm_level = line['assembly_level']
        ftp_source = line['ftp_path']

        try:
            date = datetime.date(*(int(x) for x in line['seq_rel_date'].split('/')))
        except:
            date = None

        score = asm_score(date, version_status, asm_level, False,
                          assembly_levels=self.assembly_levels)

        # store metadata for generation of assembly object later
        metadata = {
            'accession': acc,
            'bioproject': bioproj,
            'assembly_name': line['asm_name'],
            'submission_date': date,
            'version_status': version_status,
            'assembly_level': asm_level,
            'ftp_source': ftp_source,
            'tax_id': line['taxid'],
            'species_tax_id': line['species_taxid'],
            'organism_name': plain_spec,
            'infraspecific_name': plain_strain,
        }

        return bioproj, spec, strain, acc, score, rs_ref, metadata

    def read_labels(self, label_file):
        """
        read tab-seperated labelfile consisting of 2 columns containing Accession and label

        and build up a mapping index
        """
        self.mapping_index = SimpleMappingIndex()

        with open(label_file, 'r') as tab:
            for line in tab:
                line = line.strip().split('\t')
                acc = line[0]
                try:
                    label = parse_bool(line[1])
                except NoBooleanValueException as e:
                    warnings.warn(e.message)
                    continue

                p_asm = TrivialProtoAssembly(acc, label)
                self.mapping_index.add_proto_asm(p_asm)

    def read_gb_index(self):
        """
        iterate over assembly summary file (genbank) and map lines to existing organisms
        """
        valid_accessions = self.mapping_index.proto_assemblies.keys()

        with Table(self.gb_table, 'r', comment_prefix='#', skip=1) as tab:
            for i, line in enumerate(tab):
                try:
                    bioproj, spec, strain, acc, score, \
                        rs_ref, metadata = self._read_line(line, valid_accessions)

                    self.mapping_index.map(acc, score, rs_ref, metadata)
                except InvalidAssemblyException:
                    continue

    def read_rs_index(self):
        """
        iterate over assembly summary file (refseq) and map lines to existing organisms
        """
        valid_accessions = self.mapping_index.proto_assemblies.keys()

        with Table(self.rs_table, 'r', comment_prefix='#', skip=1) as tab:
            for line in tab:
                try:
                    bioproj, spec, strain, acc, score, rs_ref, \
                        metadata = self._read_line(line, valid_accessions)

                    self.mapping_index.map(acc, score, acc, metadata)
                except InvalidAssemblyException:
                    continue

    def remove_empty_organisms(self):
        """
        remove organisms that don't have at least one assembly assigned
        """
        for asm in self.mapping_index.assemblies:
            if len(asm.assemblies) == 0:
                del asm


class MappingProcessor(PreProcessor):
    """
    MappingProcessor does combined Name and Bioproject mapping. Input files have
    to be tab-seperated files, containing the columns Species, Strain, Bioproject Accession and Label

    """
    def __init__(self, label_file, mapping_log, settings, update=False):
        """
        :param str label_file: absolute path to file, containing labels
        :param bool update: update assembly summaries
        """
        self.mapping_index = None
        self.name_illegals = settings.NAME_ILLEGALS

        super(MappingProcessor, self).__init__(label_file,
                                               mapping_log, settings,
                                               update)

    def read_labels(self, label_file):
        """
        read tab-seperated labelfile consisting of columns:
        Species, Strain, Bioproject Accession, Label

        and build up a mapping index
        """

        self.mapping_index = MappingIndex()
        with Table(label_file, 'r') as tab:
            head = tab.header.header.keys()
            if not ('Species' in head and 'Strain' in head and
                    'Bioproject Accession' in head and 'Label' in head):
                print head
                raise InvalidFileFormatError
            for line in tab:
                spec, strain = unify_name(line['Species'], line['Strain'],
                                          illegals=self.name_illegals)
                bioproj = line['Bioproject Accession']
                label = None
                try:
                    label = parse_bool(line['Label'])
                except NoBooleanValueException as e:
                    warnings.warn(e.message)
                    continue

                p_asm = ProtoAssembly(bioproj, spec, strain, label)
                self.mapping_index.add_proto_asm(p_asm)

    def read_gb_index(self):
        """
        iterate over assembly summary file (genbank) and map lines to existing organisms
        """
        with Table(self.mapping_log, 'w', header=Header(
                sep='\t',
                list_=['sum.SpeciesName_plain',
                       'sum.Strain_plain',
                       'sum.SpeciesName_unif',
                       'sum.StrainName_unif',
                       'sum.Accession', 'sum.Bioproject',
                       'labelfile.SpeciesName_plain',
                       'labelfile.Strain_plain',
                       'labelfile.SpeciesName_unif',
                       'labelfile.Strain_unif',
                       'labelfile.Bioproject',
                       ])) as log:
            log.write_header()
            # TODO: rethink / check logging thing
            with Table(self.gb_table, 'r', comment_prefix='#', skip=1) as tab:
                for line in tab:
                    try:
                        bioproj, spec, strain, acc, score, \
                            rs_ref, metadata = self._read_line(line, None)
                    except InvalidAssemblyException:
                        continue

                    mapped = self.mapping_index.map(bioproj, spec, strain,
                                                    acc, score, rs_ref, metadata)

                    for asm in mapped:
                        if asm:
                            log_line = Line(log.header, [metadata['organism_name'],
                                                         metadata['infraspecific_name'],
                                                         spec, strain,
                                                         acc, bioproj, asm.plain_names[0],
                                                         asm.plain_names[1], asm.species_name,
                                                         asm.strain, asm.bioproj_acc])
                            log.write(log_line)

    def read_rs_index(self):
        """
        iterate over assembly summary file (refseq) and map lines to existing organisms
        """
        valid_accessions = self.mapping_index.rs_refs.keys()
        with Table(self.rs_table, 'r', comment_prefix='#', skip=1) as tab:
            for line in tab:
                try:
                    bioproj, spec, strain, acc, score, \
                        rs_ref, metadata = self._read_line(line, valid_accessions)
                except InvalidAssemblyException:
                    continue

                for p_asm in self.mapping_index.rs_refs[acc]:
                    p_asm.add_asm(acc, score, metadata)

    def remove_empty_organisms(self):
        """
        remove organisms that don't have at least one assembly assigned
        """
        for asm in self.mapping_index.assemblies:
            if len(asm.assemblies) == 0:
                del asm
