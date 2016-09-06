# coding: utf-8
"""
Interface to PhyloPhlAn.
"""

import logging
import os
import shutil
import subprocess
import time
from collections import deque

from phylabelle.fileio import move_files

PROJECT_DIRECTORY = os.path.split(os.getcwd())[1]


def rename_redundant_proteins(directory):
    """
    walk through all proteome-files and look for redundant protein-ids. If redundant proteins are found,
    they get renamed, by adding a suffix "_#" to their names, where "#" indicates a count. This is necessary for
    the PhyloPhlAn-pipeline to work.
    :param str directory:
    """
    print "renaming redundant proteins... ",

    id_count = {}

    files = os.listdir(directory)

    for fn in files:
        fn = os.path.join(directory, fn)
        if os.path.isdir(fn):
            continue

        str_buf = deque([])

        with open(fn, 'r+') as f:
            for line in f:
                rename = False
                if line.startswith('>'):
                    pos = line.find(' ')
                    id_ = line[1:pos]
                    try:
                        id_count[id_] += 1
                        rename = True
                    except KeyError:
                        id_count[id_] = 1

                    if rename:
                        split = line.split(' ')
                        split[0] += '_{}'.format(id_count[id_])
                        line = ' '.join(split)
                str_buf.append(line)

            f.seek(0)

            while str_buf:
                f.write(str_buf.popleft())

            f.truncate()

    print 'done'


def undo_renaming(directory):
    """
    remove the suffixes previously added by rename_redundant_proteins to restore the original IDs
    :param str directory: directory holding proteomes
    """
    print 'undo renaming... ',

    files = os.listdir(directory)

    for fn in files:
        fn = os.path.join(directory, fn)

        str_buf = deque([])

        with open(fn, 'r+') as f:
            for line in f:
                if line.startswith('>'):
                    pos1 = line.find('.')
                    pos2 = line[pos1:].find(' ') + pos1
                    if '_' in line[pos1:pos2]:
                        split_pos = line[pos1:pos2].find('_') + pos1
                        line = ''.join((line[:split_pos], line[pos2:]))
                str_buf.append(line)

            f.seek(0)

            while str_buf:
                line = str_buf.popleft()
                f.write(line)
            f.truncate()

    print 'done'


class Phylophlan(object):
    def __init__(self, settings, n_proc=2):
        """
        Init connection to PhyloPhlAn.
        :param int n_proc: number of processes, default is 2.
        """
        logging.basicConfig(filename=settings.LOGFILE, level=settings.VERBOSE_LEVEL)

        self.proc_dir = settings.PHYLOPHLAN
        self.proj_root = os.getcwd()
        self.proteome_dir = os.path.join(self.proj_root, settings.DIRECTORIES['proteome_seq'])
        self.n_proc = n_proc

        self.location = os.path.join(self.proc_dir, 'phylophlan.py')
        self.target_dir = os.path.join(self.proj_root, settings.DIRECTORIES['phylo'])
        assert os.path.isfile(self.location), 'No installation of PhyloPhlAn found'

        self.logger = logging.getLogger(__name__)

    def __call__(self, keep_tmp=False):
        """
        run complete process
        :param bool keep_tmp: keep temporary files (alignments, etc.). Note that repeated calls of PhyloPhlAn
                              without an intermediate cleanup result in an exception.
        """
        rename_redundant_proteins(self.proteome_dir)
        with self:
            self._proc()
            self._get_output()
            if not keep_tmp:
                self._cleanup()

    def __enter__(self):
        """
        Creates a symlink to use proteomes/fasta as input for PhyloPhlAn
        and changes working directory to PhyloPhlAn's location.
        """
        self.cwd = os.getcwd()

        self.symlink = os.path.join(self.proc_dir, 'input', PROJECT_DIRECTORY)

        try:
            os.symlink(self.proteome_dir, self.symlink)
        except OSError:
            pass

        assert os.path.islink(self.symlink), 'No symlink created.'

        # change to phylophlan installation directory (necessary)
        os.chdir(self.proc_dir)

    def __exit__(self, exception_type, exception_val, trace):
        """
        return to original working directory
        """
        try:
            # change back to working directory
            os.chdir(self.cwd)

            os.remove(self.symlink)
        except:
            return True

    def _proc(self):
        """
        Call PhyloPhlAn and log output.
        """
        # call phylophlan and open pipe
        proc = subprocess.Popen([self.location, PROJECT_DIRECTORY, '-u', '--nproc',
                                 str(self.n_proc)], stdout=subprocess.PIPE)

        while proc.returncode is None:
            line = proc.stdout.read()
            line.strip()
            if len(line) > 0:
                self.logger.info(line)

            if proc.stderr:
                self.logger.error(proc.stderr.read())

            proc.poll()

    def _get_output(self):
        """
        Collect phylogenetic tree.
        """
        # get output directory
        output_dir = os.path.join(self.proc_dir, 'output', PROJECT_DIRECTORY)

        # move output files to data directory
        files = os.listdir(output_dir)
        source_files = {x: os.path.join(output_dir, x) for x in files}
        # if the folder, which is supposed to contain output files already contains some older
        # version, move this to a backup folder

        # get complete paths for files in target-dir
        target_dir_files = [os.path.join(self.target_dir, x) for x in os.listdir(self.target_dir)]

        if len(target_dir_files) > 0:
            time_stamp = time.ctime(os.path.getctime(target_dir_files[0]))
            move_files(self.target_dir, 'backup_{}'.format(time_stamp))

        target_files = {x: os.path.join(self.target_dir, x) for x in files}

        for _file in files:
            shutil.move(source_files[_file], target_files[_file])

    def _cleanup(self):
        """
        calls phylophlan -c PROJECTNAME
        """
        subprocess.call([self.location, PROJECT_DIRECTORY, '-c'])

    def cleanup(self):
        """
        Delete all temporary files.
        """
        with self:
            self._cleanup()
