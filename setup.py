""" Setuptools-based setup module for bamboo

derived from the pypa example, see https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages
import os, os.path

here = os.path.abspath(os.path.dirname(__file__))

import contextlib
@contextlib.contextmanager
def chdirback(whereTo):
    import os
    cwd = os.path.abspath(os.getcwd())
    import os
    os.chdir(whereTo)
    yield
    os.chdir(cwd)

import setuptools.command.build_clib
class build_clib(setuptools.command.build_clib.build_clib):
    """ Build libraries with {"cmake" : sourcedir} in their build_info by running cmake in that directory """
    def build_libraries(self, libraries):
        for (lib_name, build_info) in libraries:
            if build_info.get("cmake"):
                self.build_cmake(lib_name, build_info.get("cmake"))
        super(build_clib, self).build_libraries([ (lib_name, build_info) for lib_name, build_info in libraries if not build_info.get("cmake") ])

    def build_cmake(self, name, cmPath):
        import pathlib
        pathlib.Path(self.build_temp).mkdir(parents=True, exist_ok=True)
        install_prefix = os.path.dirname(os.path.abspath(self.build_clib))
        buildType = "Debug" if self.debug else "Release"
        with chdirback(self.build_temp):
            self.spawn(['cmake',
                "-DCMAKE_INSTALL_PREFIX={0}".format(install_prefix),
                "-DCMAKE_BUILD_TYPE={0}".format(buildType),
                os.path.join(here, cmPath)])
            if not self.dry_run:
                self.spawn(["cmake", "--build", ".", "--target", "install", "--config", buildType])

# Get the long description from the relevant file
from io import open
with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="bamboo",

    version="0.1.0a2",

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
    scripts=[ os.path.join(root, item) for root, subFolder, files in os.walk("scripts") for item in files ],

    install_requires=["PyYAML"],
    extras_require={"Interactive mode" : ["IPython"],
                    "Tests"            : ["pytest"]},

    entry_points={},

    libraries=[("BambooExt", { "cmake" : "ext", "sources" : [ os.path.join(root, item) for root, subFolder, files in os.walk("ext") for item in files ] })],
    cmdclass={"build_clib":build_clib},
    headers=[ os.path.join(root, item) for incPath in ("cpp", os.path.join("ext", "include"))
                for root, subFolder, files in os.walk(incPath) for item in files ],
)
