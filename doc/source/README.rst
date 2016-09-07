==========
phylabelle
==========

Identify phylogenetically close pairs of proteomes, located in different
partitions across the phylogenetic tree.

What is phylabelle?
-------------------

Phylabelle is basically a pipeline, for the identification of phylogenetically close
pairs of oppositely labeled proteomes in a phylogenetic tree. It comes with an example project,
containing labels and a tree for pathogenic and non-pathogenic bacteria.

There are two purposes, one can use phylabelle for:

1.  obtain pairs from the example project, i.e. find out pairs of closely
    related pathogenic and non-pathogenic bacteria
2.  given labels, indicating some other bacterial feature, phylabelle can be used
    for the complete process: automatic download of proteomes,
    construction of phylogeny and finally obtaining of pairs.

Installation
------------

To obtain the installation files just download the the source distribution from

or just use

    ``git clone REPO``

move to phylabelle's root directory and then install it with

    ``python setup.py install``

If the installation fails because of lxml (e.g. because you're trying to install
phylabelle in a virtual environment) try it with libxml2-dev and libxslt-dev
installed.

    ``apt-get install libxml2-dev libxslt-dev``

If you intend to construct your own phylogeny, make sure, that PhyloPhlAn (see
https://huttenhower.sph.harvard.edu/phylophlan) is installed and working.
Therefore several dependencies must be satisfied:

* muscle version v3.8.31 or higher must be present in the system path and
  called "muscle"
* usearch version v5.2.32 (notice that version 6 is currently NOT supported)
  must be present in the system path and called "usearch"
* FastTree version 2.1 or higher must be present in the system path and
  called "FastTree"
* phylophlan needs some python packages that are no direct dependencies of
  phylabelle itself, and therefore not installed automatically. If you're not
  intending to construct your own phylogenies, they are not neccessary.
  Dependencies that are not covered by phylabelle are:

  - numpy
  - scipy
  - biopython


PhyloPhlAn itself is only a script and is supposed to be working without any
further installation effort. Make sure, that it is stored in a location with
writing permission.
To be used by phylabelle it is then necessary to set the entry PHYLOPHLAN in
your local_settings.py (see make your own project) to the absolute path
of the installation directory, e.g.

.. code-block:: python

    PHYLOPHLAN = '/home/USER/PATH/TO/nsegata-phylophlan-f2d78771d71d/'

Configuration
-------------

Phylabelle can be configured locally for every project. Therefore just modify
the file ``local_settings.py`` in your project folder. Make sure to use valid
python syntax. Modifiable options are explained :ref:`here<settings>`.


Operating instructions
----------------------

To init an empty project use

    ``phylabelle init_project PROJECT_NAME``

this initializes a folder structure for all datafiles generated during the
process. Note that all following phylabelle commands assume that the
project-folder is your current working directory. As a second step, use your
label-file and call

    ``phylabelle add LABEL_FILE``

LABEL_FILE must be a tab separated file, containing an assembly-accession for
refseq or genbank and a binary label, written as 0/1 or True/False (the string
gets parsed case-insensitive).
Alternatively, if the specific assemblies for each organism are unknown,
phylabelle can do a mapping-process based on organism-names and bioproject ids.
Therefore you need to provide a file, containing the fields "Species",
"Strain", "Bioproject Accession" and "Label". Then call

    ``phylabelle add -c LABEL_FILE``

to start the mapping and the download procedure.

If you don't need to construct a new phylogenetic tree, you can use the
``--no-downloads`` flag. This adds the Assemblies, without downloading
sequence files.

If you already have some proteomes stored on your device, just copy them into
the folder data/proteomes/fasta and, in case they match with the names of the
assemblies ob NCBIs ftp-site, they won't be downloaded redundantly. Attention:
Make sure, the proteome names really match the files on NCBI (have a look at
the assembly-summary files, i.e.
ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt or
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt - the
proteome files are named with [assembly_accession]_[asm_name]_protein.faa), or
leave the download to phylabelle completely. Redundant proteomes in the data
folder increase the calculation expenses for the phylogenetic tree.

Also make sure, to really call phylabelle add in order to register your
organisms and labels to phylabelles internal database.

For details on the mapping-process initiated by calling -c in phylabelle
add see the documentation under mapping.

To finally construct your phylogenetic tree use

    ``phylabelle phylophlan -n N``

where N is the number of processes you want to use for construction. Note that
a higher number of processes will speed up the construction significantly, if
your machine is able to handle it. Note that PhyloPhlAn is not reentrant.


Example project
---------------

The example project uses a metadata table from JGIs IMG
(https://img.jgi.doe.gov/cgi-bin/m/main.cgi) and a set of rules to derive
labels for pathogenicity. The automatically dervived labels have been rechecked
manually. Since the specific assemblies for each set of metadata are not
recorded in the table the lines had to be matched to assembly-ids, using
the assembly-summary files, resp. the -c flag on phylabelle add. The
input file can be derived using

``python make_labels.py data/example_data.csv \``
``"Species" "Strain" "Bioproject Accession" "Label" -s "|" \``
``-c "Label" "Label Curation" > data/labels.tsv``

The make_labels.py script can be found in the example-project folder. The same
call can be made by just calling ``./call`` (and of course
``chmod +x ./call`` in advance).
The meaning of the arguments is that the columns "Species" "Strain"
"Bioproject Accession" and "Label" should be displayed, that the
column-separator is a "|" and that, if the field "Label Curation" contains
some label, that "Label" should be overwritten with that. In the end
everything is piped into a file called labels.tsv in the folder data.

Given "labels.tsv" as an input file it is now possible for phylabelle to find
the best corresponding assemblies for each line. Therefore use

    ``phylabelle add -c data/labels.tsv``

But be careful: If you only want to try it out, you might want to limit the
number of proteomes to a reasonable amount. Else phylabelle add will
download over 3700 proteomes to your computer, and due to that the construction
of the corresponding phylogeny will probably take quite a long time, depending
on your resources of course. E.g. you could call

    ``head -n 20 data/labels.tsv > data/labels_short.tsv``

to limit the number of proteomes to 20 and then

    ``phylabelle add -c data/labels_short.tsv``

Requesting results
------------------

After having a phylogeny, phylabelle can derive pairs. Two general modes
are supported by the algorithm:

1. Derive all close pairings between members of the different groups, that
have a distance that is lower than a certain threshold value.

    ``phylabelle get_pairs -a THRESHOLD``


2. Show the closest neighbor for a specific query organism.

    ``phylabelle get_pairs -q ACCESSION``


To show all available accessions use

    ``phylabelle ls``

Worth mentioning is the ``-b`` switch of get_pairs: using this, phylabelle
will perform a minimum distance matching between all pairs in the two groups.
A plain call to ``phylabelle get_pairs`` is equivalent to
``phylabelle get_pairs -a inf -b``.
