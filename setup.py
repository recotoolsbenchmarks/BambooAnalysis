""" Setuptools-based setup module for bamboo """

from setuptools import setup, find_packages
import os, os.path
import pkg_resources

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

import setuptools.command.build_clib, setuptools.command.install
class BuildClibCommand(setuptools.command.build_clib.build_clib):
    """ Build libraries with {"cmake" : sourcedir} in their build_info by running cmake in that directory """
    def build_libraries(self, libraries):
        for (lib_name, build_info) in libraries:
            if build_info.get("cmake"):
                self.build_cmake(lib_name, build_info.get("cmake"))
        super(BuildClibCommand, self).build_libraries([ (lib_name, build_info) for lib_name, build_info in libraries if not build_info.get("cmake") ])

    def build_cmake(self, name, cmConfig):
        import pathlib
        pathlib.Path(self.build_temp).mkdir(parents=True, exist_ok=True)
        install_prefix = os.path.dirname(os.path.abspath(self.build_clib))
        buildType = "Debug" if self.debug else "Release"
        cmake_opts = [
                "-DCMAKE_INSTALL_PREFIX={0}".format(install_prefix),
                "-DCMAKE_BUILD_TYPE={0}".format(buildType),
                ]
        # pick up packages installed directly in the virtual env (e.g. Tensorflow-C), if loaded
        if "VIRTUAL_ENV" in os.environ:
            cmake_opts.append("-DCMAKE_PREFIX_PATH={0}".format(os.environ["VIRTUAL_ENV"]))
        try: # pick up libtorch, if torch is installed
            torch_cmake_locPath = os.path.join("share", "cmake", "Torch")
            if pkg_resources.resource_exists("torch", torch_cmake_locPath) and pkg_resources.get_distribution("torch").parsed_version >= pkg_resources.parse_version("1.2.0"):
                cmake_opts.append("-DTorch_DIR={0}".format(pkg_resources.resource_filename("torch", torch_cmake_locPath)))
        except ImportError as ex:
            print("Warning: optional dependency torch (>=1.2.0) not found, please install to enable inference with libtorch")
        with chdirback(self.build_temp):
            self.spawn(['cmake']+cmake_opts+[os.path.join(here, cmConfig["source_dir"])])
            if not self.dry_run:
                self.spawn(["cmake", "--build", ".", "--target", "install", "--config", buildType])

class InstallCommand(setuptools.command.install.install):
    """ Put the headers next to the python module """
    def finalize_options(self):
        ret = super(InstallCommand, self).finalize_options()
        self.install_headers = os.path.join(self.install_purelib, "bamboo", "include")
        return ret

try:
    from sphinx.setup_command import BuildDoc
except ImportError as ex:
    BuildDoc = None
    print("Warning: optional dependency sphinx not found, please install to build the documentation")

# Get the long description from the relevant file
from io import open
with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

projName = "bamboo"
projVersion = "0.1.0b2"
setup(
    name=projName,
    version=projVersion,

    description="A high-level HEP analysis library",
    long_description=long_description,
    long_description_content_type="text/markdown",

    url="https://gitlab.cern.ch/cp3-cms/bamboo",

    author="Pieter David",
    author_email="pieter.david@cern.ch",

    license="GPL-3.0-or-later",

    classifiers=[
        'Development Status :: 3 - Beta',

        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',

        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    keywords='ROOT DataFrame',

    packages=find_packages(".", exclude=["ext", "tests", "examples"]),

    setup_requires=["pytest-runner"],
    install_requires=["PyYAML", "numpy"],
    extras_require={"Interactive mode" : ["IPython"]},
    tests_require=["pytest"],

    entry_points={
        'console_scripts': [
            "bambooRun=bamboo.scripts.bambooRun:main",
            "makePUReWeightJSON=bamboo.scripts.makePUReWeightJSON:main",
            "checkBambooCMSJetDatabaseCaches=bamboo.jetdatabasecache:checkCMS_CLI"
            ]
        },

    libraries=[("BambooExt", { "cmake" : {"source_dir": "ext"}, "sources" : [ os.path.join(root, item) for root, subFolder, files in os.walk("ext") for item in files ] })],
    cmdclass={
        "build_clib" : BuildClibCommand,
        "install" : InstallCommand,
        "build_sphinx" : BuildDoc
        },
    command_options={
        "build_sphinx" : {
            "project": ("setup.py", projName),
            "version": ("setup.py", ".".join(projVersion.split(".")[:2])),
            "release": ("setup.py", projVersion),
            "source_dir": ("setup.py", "doc"),
            "build_dir" : ("setup.py", "doc/build")
            }
        },
    headers=([ os.path.join(root, item) for incPath in ("cpp", os.path.join("ext", "include"))
                for root, subFolder, files in os.walk(incPath) for item in files ]
            +[ os.path.join("build", "include", os.path.basename(fnm.strip()))
                for fnm in open(os.path.join(here, "ext", "jetclasses_filenames.txt"))
                if os.path.basename(os.path.dirname(fnm)) == "interface" ]),
)
