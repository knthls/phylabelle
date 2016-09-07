==========
phylabelle
==========

Identify phylogenetically close pairs of proteomes, located in different
partitions across the phylogenetic tree.

Documentation can be found `here<http://phylabelle.readthedocs.io/en/latest/>`_

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

To obtain the installation files just download the zip from
https://github.com/knthls/phylabelle or if you prefer using git just enter

	``git clone https://github.com/knthls/phylabelle.git``

in a terminal window. Then move to phylabelles root directory and install it with

	``python setup.py install``

If the installation fails because of lxml (e.g. because you're trying to install
phylabelle in a virtual environment) try it with libxml2-dev and libxslt-dev
installed.

    ``apt-get install libxml2-dev libxslt-dev``
	
If you intend to construct your own phylogeny, taking advadntage of the built in
wrapper for PhyloPhlAn, make sure PhyloPhlAn is installed and working. (see
https://huttenhower.sph.harvard.edu/phylophlan or check out the tutorial
`here<http://phylabelle.readthedocs.io/en/latest/phylophlan.html>`_)

In case you need to install phylabelle locally, just use

	``python setup.py install --user``

and make sure, that ``~/.local/bin`` is part of your ``$PATH``-Variable.
