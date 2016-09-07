#!/usr/bin/env python
import sys
from distutils.util import convert_path

from setuptools import setup, find_packages

if sys.version_info[:2] < (2, 7) or sys.version_info[0] > 3:
    raise RuntimeError("Python version 2.7 required.")

main_ns = {}
with open(convert_path('phylabelle/version.py'), 'r') as ver_file:
    exec(ver_file.read(), main_ns)

with open(convert_path('README.rst'), 'r') as rm_file:
    readme = rm_file.read()

setup(
    name='phylabelle',
    version=main_ns['__version__'],
    description='Identify phylogenetically close pairs of proteomes, located '
                'in different partitions across the phylogenetic tree.',
    author='Christian Knauth',
    author_email='christian.knauth@fu-berlin.de',
    packages=find_packages(),
    install_requires=['nose>=1.3.4',
                      'ete2>=2.2',
                      'sortedcontainers>=1.4.4',
                      'networkx>=1.9.1',
                      'lxml',
                      'sqlalchemy>=1.0.4',
                      'tabulate>=0.7.5',
                      ],
    entry_points={
        'console_scripts': [
            'phylabelle = phylabelle.ui:run',
        ]
    },
    include_package_data=True,
    package_data={
        'phylabelle': [
            'example_project/data/example_data.tsv',
            'example_project/make_labels.py',
            'example_project/WORKFLOW.txt',
            'example_project/sample_call'
        ],
    },
    license='GNU GPL',
    long_description=readme,
)
