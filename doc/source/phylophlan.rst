=====================
Installing PhyloPhlAn
=====================

Phylabelle is explicitly designed to use `PhyloPhlAn <https://huttenhower.sph.harvard.edu/phylophlan>`_ as a pipeline 
to generate phylogenetic models from scratch. However PhyloPhlAn is not a one-command-ready-to-install tool since 
it is based on some non-free dependencies.

Here, I will describe a way to install all dependencies without root-permissions and integrate PhyloPhlan with Phylabelle. 
This method is tested on Ubuntu 16.04 and 12.04, and should work in a similar manner on most other Linux distributions.

First we need a folder to store the executable binaries. A good folder to use for this purpose is ``~/.local/bin``, but any other 
directory will do. To tell the system, that this is going to be a folder, where it should look for executable files, we 
have to edit ``~/.bashrc`` and add the line ``PATH=$PATH:~/.local/bin``. Of course, if you're using a different directory 
the entry needs to be changed accordingly. 

Now download `muscle <http://www.drive5.com/muscle/>`_ , `usearch <http://www.drive5.com/usearch/>`_ (make sure to use v5.2.32 - 
newer versions won't work) and `FastTree <http://www.microbesonline.org/fasttree/#Install>`_,unpack and rename the executables 
to *usearch*, *FastTree* and *muscle*. Finally move the files to the folder you just created. Make sure all of them are executable with 
``chmod +x FastTree muscle usearch``.

Aditionally PhyloPhlan needs scipy, numpy and biopython. To install these, open a terminal and enter ``pip install --user numpy``, ``pip install --user scipy`` and
``pip install --user biopython``.

To install PhyloPhlAn itself just download it from `here <https://huttenhower.sph.harvard.edu/phylophlan>`_, unpack it, and edit the 
:ref:`PHYLOPHLAN <settings-phylophlan>` entry in the local settings of your project. Make sure that you have writing permissions for PhyloPhlAn's location.
