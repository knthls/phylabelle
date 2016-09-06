#!/usr/bin/env python

import gzip
import os
import sys
import warnings
from cStringIO import StringIO
from collections import defaultdict
from ftplib import FTP

from phylabelle.orm import Assembly
from phylabelle.utils import parse_bool


class Buffer(object):
    """
    Wraps a StringIO object as Buffer, to store a binary download stream and allows
    to decompress it on the fly.
    """

    def __init__(self):
        self.buf = StringIO()

    def store(self, stuff):
        self.buf.write(stuff)

    def write_and_decompress(self, path):
        self.buf.seek(0)
        decompressed = gzip.GzipFile(fileobj=self.buf, mode='rb')

        with open(path, 'wb') as f:
            f.write(decompressed.read())

    def write(self, path):
        self.buf.seek(0)

        with open(path, 'w') as f:
            f.write(self.buf.read())

    def clear(self):
        self.buf.seek(0)
        self.buf.truncate()

    def close(self):
        self.buf.close()


class FTPConnector(object):
    """
    Wraps the FTP connection from ftplib into an object, specialized for NCBIs ftp-repository.
    """

    def __init__(self, url, user, passwd):
        """
        :param str url: base directory url
        """
        self.url = url
        self._user = user
        self._passwd = passwd

        self.trace = []
        self.buf = Buffer()

    def __enter__(self):
        self.ftp = FTP(self.url)
        self._login()

        return self

    def __exit__(self, type, value, traceback):
        self.buf.close()
        try:
            self.ftp.quit()
        except:
            self.ftp.close()

    def _login(self):
        self.ftp.login(user=self._user, passwd=self._passwd)

    def remote_ls(self):
        """
        ls files in cwd
        """
        return self.ftp.nlst()

    def cwd(self, _dir):
        if _dir == '..':
            try:
                self.trace.pop()
            except:
                pass
        else:
            self.trace.append(_dir)

        self.ftp.cwd(_dir)

    def get_file(self, source_filename, target_filename, decompress=False):
        """
        retrieve a file from the ftp repository
        :param str source_filename: filename / relative path to current remote directory
        :param str target_filename: rename downloaded file
        :param bool decompress: decompress gzip compressed files
        """
        cmd = 'RETR {}'.format(source_filename)
        self.buf.clear()

        if not os.path.isfile(target_filename):
            # print 'Downloading: {}'.format(source_filename)
            self.ftp.retrbinary(cmd, self.buf.store)

            if decompress:
                self.buf.write_and_decompress(target_filename)
            else:
                self.buf.write(target_filename)


def update_indices(settings):
    """
    update assembly-indices from ncbi
    """
    print 'Updating assembly summaries... '
    rs_tab = os.path.join('data', settings.REFSEQ_TABLE)
    gb_tab = os.path.join('data', settings.GENBANK_TABLE)

    rs_source = settings.FTP_SOURCES[settings.REFSEQ_TABLE]
    gb_source = settings.FTP_SOURCES[settings.GENBANK_TABLE]

    with FTPConnector(settings.FTP_REPOSITORY, user='anonymous',
                      passwd=settings.EMAIL) as conn:
        conn.get_file(rs_source, rs_tab)
        print '    {}\t[OK]'.format(rs_tab)
        rs_tab = os.path.abspath(rs_tab)

        conn.get_file(gb_source, gb_tab)
        print '    {}\t[OK]'.format(gb_tab)
        gb_tab = os.path.abspath(gb_tab)

    return gb_tab, rs_tab


def list_dir(dir_):
    d = defaultdict(list)
    for name in os.listdir(dir_):
        # isolate assembly accession from filename
        key = '_'.join(name.split('_')[:2])
        path = os.path.join(dir_, name)
        d[key].append(os.path.join(path))
    return d


class IncompleteAssemblyException(Exception):
    def __init__(self):
        pass


class RemoteError(Exception):
    def __init__(self):
        pass


def update_labels(filename, session):
    """
    Update labels
    :param filename: file, containing new labels
    :param session: session
    :return: None
    """
    with open(filename, 'r') as f:
        for line in f:
            acc, label = line.strip().split('\t')
            label = parse_bool(label)

            asm = session.query(Assembly).get(acc)
            asm.label = label

            session.add(asm)

    session.commit()


def add_assemblies(source, session):
    """
    adds assemblies to db withou triggering any downloads
    :param source: iterable of proto assemblies
    :param session: sqlalchemy session
    """
    existing_accessions = set([x[0] for x in session.query(Assembly.accession).all()])
    for proto_asm in source:
        try:
            score, acc, metadata = proto_asm.assemblies[0]
        except IndexError:
            continue

        asm = Assembly.from_proto(proto_asm, metadata, [])
        if asm.accession not in existing_accessions:
            session.add(asm)
            existing_accessions.add(asm.accession)

    session.commit()


class Downloader(object):
    def __init__(self, source, target_dir, session, settings, force_overwrite=False):
        """
        :param list source: list of proto_assemblies with mapped assemblies
        :param bool force_overwrite: if true, all files will be downloaded, regardless whether
                                     they already exist or not, i.e. if this is set to False,
                                     the download process will be reentrant
        """
        self.source = source
        self.target_dir = target_dir
        self.session = session
        self.connector = FTPConnector(settings.FTP_REPOSITORY, user='anonymous',
                                      passwd=settings.EMAIL)
        self.existing_seqfiles = list_dir(os.path.join(target_dir, settings.DIRECTORIES['proteome_seq']))
        self.existing_metafiles = list_dir(os.path.join(target_dir, settings.DIRECTORIES['proteomes']))
        self.force_overwrite = force_overwrite
        self.settings = settings

    def run(self):

        existing_accessions = set([x[0] for x in self.session.query(Assembly.accession).all()])

        total = len(self.source)

        sys.stdout.write('Downloading proteomes\n')
        success = 0
        omitted = 0

        with self.connector as conn:
            for i, proto_asm in enumerate(self.source):
                try:
                    for score, acc, metadata in proto_asm.assemblies:

                        try:
                            # get remote directory
                            dir = metadata['ftp_source'].split('.gov')[1]
                        except IndexError:
                            # try next asm
                            continue

                        if not (any(self.check_components(acc)) or self.force_overwrite):
                            if acc not in existing_accessions:
                                paths = self.existing_seqfiles[acc] + self.existing_metafiles[acc]
                                asm = Assembly.from_proto(proto_asm, metadata, paths)
                                self.session.add(asm)
                                existing_accessions.add(asm.accession)

                            break

                        try:
                            # change into assembly directory on remote
                            conn.cwd(dir)
                        except EOFError:
                            raise RemoteError

                        # list available files
                        files = conn.remote_ls()
                        # choose files which are supposed to be downloaded
                        files = [x for x in files if ('assembly_report' in x) or
                                 ('protein.faa' in x) or ('feature' in x)]

                        if len(files) < 3:
                            continue

                        try:
                            paths = []

                            for f in files:
                                if 'protein.faa' in f:
                                    target_dir = os.path.join(
                                        self.target_dir,
                                        self.settings.DIRECTORIES['proteome_seq']
                                    )
                                else:
                                    target_dir = os.path.join(
                                        self.target_dir,
                                        self.settings.DIRECTORIES['proteomes']
                                    )

                                if '.gz' in f:
                                    target_filename = os.path.join(target_dir, f[:-3])
                                    decompress = True
                                else:
                                    target_filename = os.path.join(target_dir, f)
                                    decompress = False

                                if os.path.isfile(target_filename) and not self.force_overwrite:
                                    paths.append(target_filename)
                                    continue

                                try:
                                    conn.get_file(f, target_filename, decompress=decompress)

                                    paths.append(target_filename)
                                except Exception as e:
                                    warnings.warn(e)
                                    raise IncompleteAssemblyException

                            asm = Assembly.from_proto(proto_asm, metadata, paths)
                            if asm.accession not in existing_accessions:
                                self.session.add(asm)
                                existing_accessions.add(asm.accession)

                            break

                        except IncompleteAssemblyException:
                            continue

                    success += 1
                except RemoteError:
                    omitted += 1
                    sys.stderr.write('Remote server error. Ending download. Please restart procedure manually.')
                    break
                except KeyboardInterrupt as e:
                    omitted += 1
                    sys.stderr.write('Interrupting download procedure; Progress gets committed')
                    break

                sys.stdout.write('    {} of {}, found: {}, omitted: {}'.format(i+1, total, success, omitted))
                if i < total - 1:
                    sys.stdout.write('\r')
                    sys.stdout.flush()
                else:
                    sys.stdout.write('\n')

        self.session.commit()

    def check_components(self, key):
        """
        for a given Assembly-id, check whether sequence-file, assembly-report and feature-table exist.
        returns 3 booleans, one for each, where False indicates, that a file exists, True
        that it's missing
        :param str key: assembly accession
        """
        seqfile = len(self.existing_seqfiles[key]) == 0
        report = not any(['report' in x for x in self.existing_metafiles[key]])
        feature_table = not any(['feature' in x for x in self.existing_metafiles[key]])

        return seqfile, report, feature_table
