#!/usr/bin/env python2.7

from setuptools import setup, find_packages

setup(
    name='EMBLmyGFF3',
    version='1.2.2',

    description='An efficient way to convert gff3 annotation files into EMBL format ready to submit',

    url='https://github.com/NBISweden/EMBLmyGFF3',
    download_url='https://github.com/NBISweden/EMBLmyGFF3/archive/v1.2.2.tar.gz',
    author='Martin Norling, Niclas Jareborg, Jacques Dainat',

    license='GPL-3.0',
    packages=find_packages(),

    install_requires=['biopython==1.67', 'bcbio-gff==0.6.4'],
    include_package_data=True,

    entry_points={
        'console_scripts': ['EMBLmyGFF3 = EMBLmyGFF3:main',
        ],
    }
)
