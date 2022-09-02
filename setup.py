#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# access the version wihtout importing the EMBLmyGFF3 package
with open('EMBLmyGFF3/version.py') as f: exec(f.read())

setup(
    name='EMBLmyGFF3',
    version=__version__,

    description='An efficient way to convert gff3 annotation files into EMBL format ready to submit',

    url='https://github.com/NBISweden/EMBLmyGFF3',
    download_url='https://github.com/NBISweden/EMBLmyGFF3/archive/v' + __version__ +'.tar.gz',
    author='Martin Norling, Niclas Jareborg, Jacques Dainat',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['biopython>=1.78', 'bcbio-gff>=0.6.4','numpy>=1.22', 'python_version>="3.6.0"' ],
    include_package_data=True,

    entry_points={
        'console_scripts': ['EMBLmyGFF3 = EMBLmyGFF3:main',
        'EMBLmyGFF3-augustus-example = examples.augustus_example:main',
        'EMBLmyGFF3-maker-example = examples.maker_example:main',
        'EMBLmyGFF3-prokka-example = examples.prokka_example:main',
        ],
    }
)
