""" Setuptools-based setup module for bamboo

derived from the pypa example, see https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
from io import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="bamboo",

    version="0.1.0a1",

    description="A high-level HEP analysis library for ROOT::RDataFrame",
    long_description=long_description,
    long_description_content_type="text/markdown",

    url="https://cp3-git.irmp.ucl.ac.be/pdavid/bamboo",

    author="Pieter David",
    author_email="pieter.david@uclouvain.be",

    license="GPL-3.0-or-later",

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',

        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='ROOT DataFrame',

    packages=["bamboo"],

    install_requires=[],

    extras_require={},

    package_data={},
    data_files=[],

    entry_points={},
)
