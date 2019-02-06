#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, calour development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup


classifiers = [
    'Development Status :: 2 - Pre-Alpha',
    'License :: OSI Approved :: BSD License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'phenotype-calour is a calour interface to the IJSEM_pheno_db database (see "Hiding in Plain Sight: Mining Bacterial Species Records for Phenotypic Trait Information" - http://msphere.asm.org/content/2/4/e00237-17)'

with open('README.md') as f:
    long_description = f.read()

keywords = 'microbiome calour phenotype database analysis bioinformatics',

setup(name='pheno-calour',
      version=0.1,
      license='BSD',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="calour development team",
      maintainer="calour development team",
      url='https://github.com/amnona/calour',
      test_suite='nose.collector',
      packages=find_packages(),
      package_data={'phenocalour': ['data/pheno-calour-data.txt']},
      install_requires=[
          'calour'],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'coverage': ["coverage"],
                      'doc': ["Sphinx >= 1.4"]},
      entry_points={
          'console_scripts': [
              'calour=calour.cli:cmd',
          ]})
