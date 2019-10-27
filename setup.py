import sys
from distutils.core import setup, Extension
import os
import warnings

PYPINAME = "fastq-and-furious"
PACKAGENAME = "fastqandfurious"
VERSION="0.3.0"

extra_compile_args = ['-Wall']

CLASSIFIERS = [
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: C++",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering",
]

if tuple(sys.version_info[:2]) < (3, 5):
    print("Error: Python >= 3.5 is *required*.")
    sys.exit(1)
    

if sys.platform == 'darwin':
    warnings.warn("Not tested on OSX. Feedback welcome.")
    pass
elif sys.platform == 'linux':
    pass
else:
    raise ValueError("The platform %s is not supported." % sys.platform)

faf_mod = Extension("%s._fastqandfurious" % PACKAGENAME,
                    sources=["src/_fastqandfurious.c", ],
                    #depends=["src/.h"],
                    #include_dirs=["src",],
                    language="c",
                    extra_compile_args = extra_compile_args + \
                    ['-O3', '-std=c99'])

setup(
    name = PYPINAME,
    version = VERSION,
    description = "Fast handling of FASTQ files",
    license = "MIT",
    author = "Laurent Gautier",
    author_email = "lgautier@gmail.com",
    url = "https://github.com/lgautier/fastq-and-furious",
    packages = [PACKAGENAME,
                #PACKAGENAME + '.tests',
                PACKAGENAME + '.demo'],
    package_dir = {PACKAGENAME: 'src'},
    ext_modules = [faf_mod, ],
    extras_require = {
        'test' : ['pytest', ],
        'demo' : ['ngs_plumbing', 'screed', 'biopython']},
    classifiers = CLASSIFIERS)



