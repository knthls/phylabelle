"""
phylabelle.fileio
=================

Input- and output-functions
"""

import os
import shutil
import sys
import warnings

import tabulate


class AbsolutePathException(Exception):
    """
    Exception which is raised, if an absolute path instead of a relative one
    is found
    """
    def __int__(self):
        pass


def split_dirs(dir_name):
    """
    recursively iterate over a path, and return a list of directories from
    general to specific
    :param str dir_name: path
    :yield: paths
    """
    top, name = os.path.split(dir_name)
    if top == '.':
        return [dir_name]
    elif top == '/':
        raise AbsolutePathException
    else:
        return split_dirs(top) + [dir_name]


def create_project_paths(base_dir, subdirs):
    """
    Creates a project folder-structure
    :param str base_dir: directory to use
    :param subdirs: structure as configured in settings
    """

    for dir_ in subdirs.values():
        dir_ = os.path.join(base_dir, dir_)
        try:
            os.makedirs(dir_)
        except OSError:
            pass


def is_valid_project(base_dir, subdirs):
    """
    Creates a project folder-structure
    :param str base_dir: directory to use
    :param subdirs: structure as configured in settings
    """
    for dir_ in subdirs.values():
        dir_ = os.path.join(base_dir, dir_)
        if not os.path.isdir(dir_):
            return False

    return True


def move_files(source, target):
    """
    recursively move files from source directory to target directory
    :param str source: source directory
    :param str target: target directory
    """
    # make directory
    try:
        os.makedirs(target)
    except OSError:
        pass
    shutil.move(source, target)


def format_get_lines(str_, std_len=80):
    """
    Generator to reformat strings in order to get lines for a given maximum length. Removes all
    previous formatting, newlines, tabs and whitespace. Joins with \' \' .

    :param str str_: the text which should be formatted
    :param int std_len: the maximum length of the line
    :yield: formatted lines of len(line) <= std_len
    """
    c = 0
    line = []
    for word in str_.split():
        word = word.strip()
        l = len(word)
        if c + l >= std_len:
            yield ' '.join(line)
            line = []
            c = 0

        line.append(word)
        c += l

    yield ' '.join(line)


def read_phyloxml(input_file):
    """
    Parses a pyhlogenetic tree in phyloxml-format.
    :param str input_file: path to file
    :return: ete2.Tree object
    """
    from ete2 import Phyloxml
    project = Phyloxml()
    project.build_from_file(input_file)

    return project.get_phylogeny()


def read_tree(tree_file):
    """
    Read a phylogenetic tree and return as ete2.Tree object. Can handle
    newick and phyloxml.
    :param str tree_file:
    :return: ete2.Tree object
    """
    from ete2 import Tree
    filename = os.path.split(tree_file)[1]
    extension = os.path.splitext(filename)[1]

    if extension == '.xml':
        return read_phyloxml(os.path.abspath(tree_file))[0]
    else:
        try:
            return Tree(os.path.abspath(tree_file))
        except:
            raise ValueError('filetype has to be either \'nwk\' or \'xml\'')


class Header(object):
    """
    Header assumes, that a table's header line has the maximum length of any line
    """

    def __init__(self, sep, line=None, list_=None, prefix=None):
        """
        :param str sep: seperator of fields
        :param str line: string describing fields of a table, i.e. a header string
        :param list list_: provide the list fields directly (instead of line and sep)
        :param str prefix: indicates that the header string contains a prefix that \
                    does not belong to the description, e.g. "#" or other \
                    commenting characters
                    """
        assert line or list_, 'invalid initialization of Header class'

        if line:
            self.str_ = line
            if prefix:
                if line.startswith(prefix):
                    line = line[line.find(prefix) + len(prefix):]

            list_ = line.strip().split(sep)
        else:
            self.str_ = sep.join(list_)

        self.list_ = list_
        self.header = {b.strip(): a for a, b in enumerate(list_)}

    def __len__(self):
        return len(self.header)

    def __iter__(self):
        return self.header.iterkeys()

    def get_value(self, key, line):
        return line[self.header[key]]

    def get_index(self, key):
        return self.header[key]


class Line(object):
    def __init__(self, header, values=None):
        self.header = header
        if values is None:
            self.values = ['NA'] * len(header)
        else:
            self.values = values

    def __getitem__(self, key):
        return self.header.get_value(key, self.values)

    def __setitem__(self, key, value):
        self.values[self.header.get_index(key)] = value

    def __iter__(self):
        return iter(self.values)

    def consolidate(self, sep, fill=False):
        if not fill:
            assert len(self.values) == len(self.header), 'invalid line length'
            return sep.join([str(x) for x in self.values]) + '\n'
        else:
            return sep.join([str(x) for x in self.values] + [''] * (len(self.header) - len(self.values))) + '\n'


class Table(object):
    """
    Table uses the header strings of a csv-formatted table for accessing the
    value of a column in a specific Line (see `Line`). To read in the values
    of column ``key`` da something like this:

    .. code-block:: python

        with Table('filename', 'r') as tab:
            for line in tab:
                value = line['key']

    Table also supports writing-access:

    .. code-block:: python

        with Table('filename', 'w', header=Header(
                                                sep = '\t',
                                                list_=['col1', 'col2']
                                                ) as tab:
            for val1, val2 in values:
                tab.write([val1, val2])
    """
    def __init__(self, filename, access_mode,
                 header=None, seperator='\t',
                 comment_prefix=None, skip=0):
        """
        :param str filename:
        :param access_mode: can be either 'w' or 'r'
        :param Header header: `Header`-object
        :param str seperator:
        :param str comment_prefix: in case the header line begins with a
        prefix that is not part of the first column name itself, this can
        be specified here
        :param int skip: number of lines to be skipped
        """
        self.filename = filename
        self.access_mode = access_mode
        self.seperator = seperator
        self.comment_prefix = comment_prefix
        self.header = header
        self.skip = skip

    def __iter__(self):
        assert hasattr(self, 'file'), 'Needs to be entered via with statement'
        return self.read()

    def __enter__(self):
        if self.filename == 'stdout':
            self.file = sys.stdout
        else:
            self.file = open(self.filename, self.access_mode)

        if self.access_mode == 'r' and self.header is None:
            for i in xrange(self.skip):
                self.file.next()
            self.header = Header(self.seperator, line=self.file.next(),
                                 prefix=self.comment_prefix)
        else:
            if self.header is None:
                warnings.warn('no Header found')

        return self

    def __exit__(self, exception_type, exception_val, trace):
        if not self.file == sys.stdout:
            try:
                self.file.close()
            except:
                return True

    def write_header(self):
        """
        write header line to file.
        :return:
        """
        self.file.write(self.header.str_)
        self.file.write('\n')

    def read(self):
        for line in self.file:
            line = line.replace('\n', '').split(self.seperator)
            yield Line(self.header, line)

    def write(self, line):
        if type(line) == str:
            str_ = line
        else:
            str_ = self.seperator.join(line) + '\n'

        self.file.write(str_)


class BufferedTable(Table):
    """
    inherits from table, but enables to store values in a buffer
    """
    def __init__(self, filename, access_mode,
                 header=None, seperator='\t',
                 comment_prefix=None):
        super(BufferedTable, self).__init__(filename, access_mode,
                                            header, seperator, comment_prefix)

        self._buf = []

    def __enter__(self):
        assert self.filename is not None, 'unknown filename'
        assert self.access_mode is not None, 'unknown filename'
        return super(BufferedTable, self).__enter__()

    def __iter__(self):
        return iter(self._buf)

    def __getitem__(self, key):
        assert type(key) == int
        return self._buf[key]

    def append(self, line):
        # assert isinstance(line, Line)
        self._buf.append(line)

    def write(self, sort_by=None, pretty=False, consolidate=False):
        """
        :param str sort_by: Field in Header which provides a key for sorting
        :param bool pretty: boolean indicating whether to use tabulate for formatting
        :param bool consolidate: boolean indicating whether values are list or strings, \
                                i.e. if a line has to be joined to a string or not
        """
        if sort_by:
            num = self.header.header[sort_by]
            items = sorted(self._buf, key=lambda x: x[num])
        else:
            items = self._buf

        if pretty:
            self.file.write(tabulate.tabulate(items, headers=self.header.list_))
            self.file.write('\n')
        else:
            for item in items:
                if consolidate:
                    item = self.seperator.join([str(x) for x in item])
                self.file.write(item)
                self.file.write('\n')


def fmt_pair(pair):
    """
    :param pair: tuple of the form (Assembly1, Assembly2, distance)
    :return:
    """
    asm1, asm2, dist = pair
    if asm1.label:
        t_asm = asm1
        f_asm = asm2
    else:
        t_asm = asm2
        f_asm = asm1

    return [t_asm.accession, t_asm.organism_name, f_asm.accession, f_asm.organism_name, dist]


def get_pretty_output(results, sort_by):
    """
    takes a list of pairs and prints it in a pretty way
    :param results: iterable of pairs, i.e. tuples of the form *(accession1, accession2, distance)*
    :param str sort_by:
    """
    header = ['Positive Accession', 'Positive Name', 'Negative Accession', 'Negative Name', 'Distance']

    # the tuples are sorted, so set should work
    # results = set(results)

    if sort_by == 'dist':
        func = lambda line: line[4]
    elif sort_by == 'p_name':
        func = lambda line: line[1]
    else:
        func = lambda line: line[3]

    print tabulate.tabulate(sorted((fmt_pair(x) for x in results), key=func), headers=header)
