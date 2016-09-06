#!/usr/bin/env python
# coding: utf8
"""
phylabelle
==========

Identify phylogenetically close pairs of proteomes, located in different
partitions across the phylogenetic tree.

"""

import argparse
import os
import shutil
import sys
import warnings

import tabulate
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from phylabelle.core import PhyloTree, NoResultsException
from phylabelle.fileio import get_pretty_output, BufferedTable, Header, \
    create_project_paths, is_valid_project
from phylabelle.maintenance import update_labels, add_assemblies

# this enables to use local settings
sys.path.append(os.getcwd())
try:
    import local_settings as settings
except ImportError as e:
    from phylabelle import settings
sys.path.remove(os.getcwd())


def db_connect(subdir=None, init=False):
    """
    :param subdir: if sqlite-db is stored in a directory
                                  different from current working directory
    :param bool init: if true, database file will be initialized
    """
    if subdir is None:
        conn = 'sqlite:///{}'.format(settings.DB)
    else:
        conn = 'sqlite:///{}'.format(os.path.join(subdir, settings.DB))

    engine = create_engine(conn, echo=False)
    session = sessionmaker(bind=engine)()

    if init:
        init_orm(engine, session)

    return session


def copytree(src, dst):
    """
    recursively copy a folder structure (wrapper for shutil.copytree)
    :param src: source directory
    :param dst: target directory
    """
    for item in os.listdir(src):
        src_ = os.path.join(src, item)
        dst_ = os.path.join(dst, item)
        if os.path.isdir(src_):
            try:
                shutil.copytree(src_, dst_)
            except OSError:
                copytree(src_, dst_)
        else:
            shutil.copy2(src_, dst_)


def init_blank_project(directory):
    """
    create folder structure and database file
    :param str directory: subdirectory
    """
    create_project_paths(directory, settings.DIRECTORIES)
    db_connect(subdir=directory, init=True)


def init_orm(engine, session):
    """
    create database file and models
    :param engine: sqlalchemy db-engine
    :param session: sqlalchemy session
    """
    from phylabelle import orm

    orm.Base.metadata.create_all(engine)
    session.commit()
    session.close()


def print_ls_asm(assemblies, sort_by='Name', pretty=False, show_plasmids=True):
    """
    print a plain table of assemblies
    :param iter assemblies: iterable of assemblies
    :param sort_by:
    :param pretty:
    :param show_plasmids:
    :return:
    """
    # TODO: show plasmids
    header = Header('\t', list_=settings.LS_ASM_HEADER.keys())

    fmt = []

    for name, attr in settings.LS_ASM_HEADER.iteritems():
        try:
            func = eval(settings.LS_ASM_HEADER[name])
        except:
            func = None

        fmt.append((attr, func))

    with BufferedTable('stdout', 'w', header=header) as tab:
        if not pretty:
            tab.write_header()

        for asm in assemblies:
            line = []
            for attr, func in fmt:
                val = getattr(asm, attr)
                if func:
                    val = func(val)
                line.append(val)

            tab.append(line)

        tab.write(sort_by=sort_by, pretty=pretty, consolidate=True)


def print_get_pairs(results, sort_by):
    header = Header('\t',
                    list_=['Positive Accession',
                           'Positive Name', 'Negative Accession',
                           'Negative Name', 'Distance'])

    results = set(results)

    if sort_by == 'dist':
        func = lambda line: line[4]
    elif sort_by == 'p_name':
        func = lambda line: line[1]
    else:
        func = lambda line: line[3]

    print tabulate.tabulate(sorted((x for x in results), key=func), headers=header)


def print_query_pairs(query, partners):
    """
    print pairs
    """

    tup = partners[0]
    if tup[0].accession == query:
        # get assembly object
        q_asm = tup[0]
    else:
        q_asm = tup[1]

    print "query:\n======"

    print "  organism name: {:20}\n" \
          "  species name: {:20}\n" \
          "  infraspecific name: {:20}\n" \
          "  assembly accession: {:20}".format(q_asm.organism_name, q_asm.species_name,
                                              q_asm.infraspecific_name,
                                              q_asm.accession)

    print "\n"
    print "partners:\n=========\n"

    lines = []

    for p1, p2, dist in partners:
        # p1, p2 are assembly objects
        if p1.accession == query:
            line = [p2.organism_name, p2.species_name, p2.infraspecific_name, p2.accession, dist]
        else:
            line = [p1.organism_name, p1.species_name, p1.infraspecific_name, p1.accession, dist]

        lines.append(line)

    print tabulate.tabulate(lines, headers=[
        'organism name',
        'species name',
        'infraspecific name',
        'accession',
        'distance'
    ])


class FriendlyArgumentParser(argparse.ArgumentParser):
    """
    extends argparse.ArgumentParser to print the help message if no
    arguments are given.
    """

    def error(self, message):
        """
        overrides default error method to print help-message
        :param message:
        :return:
        """
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def get_args():
    """
    contains all argument parser settings
    :return: parsed arguments
    """
    #: ArgumentParser-instance
    arg_parser = FriendlyArgumentParser(description=__doc__)
    subparsers = arg_parser.add_subparsers(dest='subparser')

    #: store all subparsers as dict entries, keyed with subprogram name
    sp_dict = {}

    sp_dict['init_project'] = subparsers.add_parser('init_project',
                                                    help='init new phylabelle project')
    sp_dict['init_project'].add_argument('name', metavar='NAME', type=str,
                                         help='Project name')

    sp_dict['ls'] = subparsers.add_parser('ls', help='list available assemblies')

    # set up arguments as mutual exclusive groups
    mutex_groups = {}

    sp_dict['ls'].add_argument('-q', '--query', type=str, nargs=2,
                               metavar=('FIELD', 'QUERY'),
                               help='list all available organisms, \
                                        containing QUERY in FIELD')
    sp_dict['ls'].add_argument('--show_fields', action='store_true',
                               default=False,
                               help='show a list of valid fields for --query')

    sp_dict['ls'].add_argument('--tsv', default=False, action='store_true',
                               help='set this this flag to print the output as \
                                    tab-separated-values')

    sp_dict['ls'].add_argument('--labels', action='store_true', default=False,
                               help='prints a list of assembly accessions and their labels')

    sp_dict['ls'].add_argument('-s', '--sort_by', type=str,
                               choices=['Name', 'Accession',
                                        'Assembly Level', 'Label'],
                               default='Name',
                               help='sort output')

    sp_dict['get_example'] = subparsers.add_parser('get_example',
                                                   help=get_example.__doc__)

    sp_dict['get_example'].add_argument('-p', '--path', default=None, type=str,
                                        help='copy files to subdirectory, if specified.')

    sp_dict['add'] = subparsers.add_parser('add', help='download data from tab separated \
                                                       input file, mapping labels to \
                                                       assembly accessions')

    sp_dict['add'].add_argument('file', metavar='FILE', type=str, help='input-file')
    sp_dict['add'].add_argument('--no-download', action='store_true', default=False,
                                help='suppresses downloads.')

    mutex_groups['add'] = sp_dict['add'].add_mutually_exclusive_group()

    mutex_groups['add'].add_argument('-c', '--complex', action='store_true', default=False,
                                     help='Use a mapping procedure, to derive the best \
                                     available assembly from bioproject-accessions and \
                                     organism names')
    mutex_groups['add'].add_argument('-u', '--update', action='store_true', default=False,
                                     help='update labels only. No downloads will be triggered.')

    sp_dict['phylophlan'] = subparsers.add_parser('phylophlan', help='generate phylogenetic tree')

    sp_dict['phylophlan'].add_argument('-n', '--n_proc', type=int, default=2,
                                       help='number of processes to use for phylogeny-construction')
    sp_dict['phylophlan'].add_argument('-c', '--cleanup', action='store_true',
                                       help='delete additional files generated by \
                                            PhyloPhlAn (Alignments etc.)')
    sp_dict['phylophlan'].add_argument('-k', '--keep', action='store_true',
                                       help='keep all files generated in the phylogeny \
                                            construction process (Alignments, etc.). \
                                            Note that they prevent PhyloPhlAn from being run \
                                            again. For a repeated run, call phylabelle \
                                            phylophlan --cleanup in between')

    sp_dict['get'] = subparsers.add_parser('get_pairs', help='find pairs of oppositely \
                                                              labeled assemblies')

    sp_dict['get'].add_argument('--tsv', default=False, action='store_true',
                                help='set this this flag to print the output as \
                                    tab-separated-values')

    sp_dict['get'].add_argument('-b', default=False, action='store_true',
                                help='Show bijective relations only, meaning \
                                     that for two Assemblies A and B, \
                                     a pairing is only considered valid if \
                                     the "closest-partner"-relation (CP) \
                                     holds in both directions, i.e. \
                                     CP(A) = B <=> CP(B) = A')

    sp_dict['get'].add_argument('-s', '--sort_by', type=str,
                                choices=['p_name', 'np_name', 'dist'],
                                default='dist',
                                help='sort output by the name of the \
                                     pathogenic partner, the name of \
                                     the non-pathogenic partner or \
                                     their distance')

    sp_dict['get'].add_argument('-m', '--mode', type=str, default='all',
                                choices=['intra', 'inter'],
                                metavar='MODE',
                                help='Choose between "inter"(-species) \
                                      and "intra"(-species). Sets a \
                                      constraint on the type of pairs, \
                                      which can be found. If "intra" is \
                                      chosen, only pairs within the same \
                                      species are considered valid, and \
                                      according to this if "inter" is chosen \
                                      only pairs comprising different species \
                                      are shown.')
    sp_dict['get'].add_argument('--max', type=int, default=None,
                                metavar='N',
                                help='restricts the number \
                                     of pairs to be displayed to the N \
                                     closest.')
    mutex_groups['get_pairs'] = sp_dict['get'].add_mutually_exclusive_group()
    mutex_groups['get_pairs'].add_argument('-a', '--all', nargs='?',
                                           metavar='THRESHOLD', type=float, const=float('inf'),
                                           help='find all available pairs, with distance \
                                                smaller than THRESHOLD')

    mutex_groups['get_pairs'].add_argument('-q', '--query', type=str, nargs=1, metavar='QUERY',
                                           default=None, help='find the closest partner for a given \
                                                              QUERY (assembly accession)')
    mutex_groups['get_pairs'].add_argument('--mappings', action='store_true', default=False,
                                           help='get the annotation for the tree')

    return arg_parser.parse_args()


# ###########
# subprograms
# ###########


def ls(args):
    import phylabelle.orm as objects

    session = db_connect()

    if args.labels:
        data = session.query(objects.Assembly.accession, objects.Assembly.label).all()
        print 'accession\tlabel'

        for acc, label in data:
            print '\t'.join([str(acc), str(label)])
        return

    if args.show_fields:
        objects.Assembly.show_fields()
        return

    sort_by = lambda x: getattr(x, settings.LS_ASM_HEADER[args.sort_by])

    if args.query is not None:
        target = getattr(objects.Assembly, settings.LS_ASM_HEADER[args.query[0]])
        assemblies = sorted(session.query(objects.Assembly).filter(
            target.ilike('%{}%'.format(args.query[1]))).all(), key=sort_by)

        if len(assemblies) == 0:
            print 'No assemblies found, matching the query "{}" in field "{}".'.format(args.query[1], args.query[0])
            return
    else:
        assemblies = sorted(session.query(objects.Assembly).all(),
                            key=sort_by)

    print_ls_asm(assemblies, sort_by=args.sort_by, pretty=not args.tsv)


def init_project(args):
    """
    Init a project folder
    :param args: command line arguments, has to contain name
    """
    name = args.name

    if not name.startswith('./'):
        name = './' + name

    init_blank_project(name)

    dest = os.path.abspath(os.getcwd())
    dest = os.path.join(dest, args.name)

    shutil.copy(os.path.join(os.path.split(__file__)[0], 'settings.py'),
                os.path.join(dest, 'local_settings.py'))
    print 'project {} created'.format(name)


def get_pairs(args, phylo_tree):
    max_ = args.max

    if args.all:
        threshold = args.all

        phylo_tree.evaluate_all_pairs(threshold,
                                      mode=args.mode)
        if args.b:
            results = phylo_tree.get_minimum_matching()
        else:
            try:
                results = phylo_tree.get_closest()
            except NoResultsException:
                print 'No Pairs found. Maybe try a higher threshold.'
                return

        if max_ is not None:
            get_pretty_output(results[:max_], args.sort_by)
        else:
            get_pretty_output(results, args.sort_by)

    elif args.mappings:
        print 'AssemblyAccession\tLabel\tSpeciesTaxID'
        for asm in phylo_tree.assemblies.itervalues():
            print '\t'.join((asm.accession, str(asm.label), str(asm.species_tax_id)))
    elif args.query:
        query = args.query[0]
        if max_ is None:
            max_ = 1

        results = phylo_tree.find_closest_partner(query, mode=args.mode,
                                                  n_results=max_)

        print_query_pairs(query, results)
    else:
        warnings.warn('If no further options are supplied, only a minimum '
                      'matching of all pairs will be shown.'
                      'This is equivalent to "--all inf -b".')
        phylo_tree.evaluate_all_pairs(float('inf'),
                                      mode=args.mode)
        results = phylo_tree.get_minimum_matching()

        get_pretty_output(results, args.sort_by)


def add(args):
    """
    Start data download with tab seperated file, consisting of assembly accessions and labels
    :param args:
    """
    from phylabelle.selection import MappingProcessor, PreProcessor, InvalidFileFormatError
    from phylabelle.maintenance import Downloader

    if args.complex:
        try:
            processor = MappingProcessor(args.file, mapping_log='./mapping_log.tsv',
                                         settings=settings, update=True)
        except InvalidFileFormatError:
            print 'Invalid file format. Make sure that your file has columns ' \
                  '"Species", "Strain", "Bioproject Accession" and "Label"'
            return
        processor.read_gb_index()
        processor.mapping_index.check_integrity()
        processor.read_rs_index()
        processor.remove_empty_organisms()
    elif args.update:
        session = db_connect()
        update_labels(args.file, session)
        return

    else:
        processor = PreProcessor(args.file, mapping_log='./mapping_log.tsv',
                                 settings=settings, update=True)
        processor.read_gb_index()
        processor.read_rs_index()
        processor.remove_empty_organisms()

    # create_project_paths('.')
    assert is_valid_project(os.getcwd(), settings.DIRECTORIES), 'Setup project structure first ' \
                                                                '(call phylophlan init)'
    session = db_connect()

    if args.no_download:
        add_assemblies(processor.mapping_index.assemblies, session)
    else:
        d = Downloader(processor.mapping_index.assemblies, '.',
                       session, settings)
        d.run()


def get_example(args):
    """
    Copy an example project to working directory
    """

    print 'init example project'

    # get project base directory
    base = os.path.split(os.path.abspath(__file__))[0]

    source = os.path.join(base, 'example_project')

    assert os.path.isdir(source)

    if args.path is not None:
        dest = args.path
    else:
        dest = os.getcwd()

    init_blank_project(os.path.join(os.path.relpath(dest), 'example_project'))

    dest = os.path.abspath(dest)
    dest = os.path.join(dest, 'example_project/')

    copytree(source, dest)

    shutil.copy(os.path.join(os.path.split(__file__)[0], 'settings.py'), os.path.join(dest, 'local_settings.py'))


def phylophlan(args):
    """
    run phylophlan on all available proteomes
    """
    from phylabelle.phylo import Phylophlan

    proc = Phylophlan(n_proc=args.n_proc, settings=settings)

    if args.cleanup:
        proc.cleanup()
        return

    proc(keep_tmp=args.keep)


def run():
    """
    Main routine
    """
    # parse arguments
    args = get_args()

    # subparser contains name of subprogram
    # subprogram has to be listed in globals()
    func = globals()[args.subparser]

    if args.subparser == 'get_pairs':
        # in this case a phylogenetic tree gets loaded
        # list comprehension finds all xml-files in phylo directory
        tree_files = [os.path.join(settings.DIRECTORIES['phylo'], x)
                      for x in os.listdir(settings.DIRECTORIES['phylo'])
                      if '.xml' in x]

        assert len(tree_files) > 0, "No phylogenetic tree found"
        assert len(tree_files) < 2, "Multiple phylogenetic trees found"

        session = db_connect()

        phylo_tree = PhyloTree(tree_files[0], session=session)
        func(args, phylo_tree)
    else:
        func(args)