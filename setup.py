#!/usr/bin/env python3
"""
Setup script for EMBLmyGFF3. This tool is developed to convert and reformat the
common GFF3 format into the more stringent EMBL format used for data submission
to the European Nucleotide Archive (ENA).
"""

from setuptools import setup, find_packages

setup(
    name='EMBLmyGFF3',
    version='2.0.0',

    description=('An efficient way to convert gff3 annotation files into EMBL '
                 'format ready to submission'),

    url='https://github.com/NBISweden/EMBLmyGFF3',
    download_url='https://github.com/NBISweden/EMBLmyGFF3/archive/v2.0.0.tar.gz',
    author='Martin Norling, Niclas Jareborg, Jacques Dainat',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['biopython==1.67', 'bcbio-gff==0.6.4'],
    include_package_data=True,

    entry_points={
        'console_scripts': [
            'EMBLmyGFF3 = EMBLmyGFF3.convert:main',
        ],
    }
)
