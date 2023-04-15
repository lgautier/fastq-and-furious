import sys
from setuptools import setup, Extension
import os
import warnings

PACKAGENAME = 'fastqandfurious'

extra_compile_args = ['-Wall']

if tuple(sys.version_info[:2]) < (3, 7):
    print('Error: Python >= 3.7 is *required*.')
    sys.exit(1)
    

if sys.platform == 'darwin':
    warnings.warn('Not tested on OSX. Feedback welcome.')
    pass
elif sys.platform == 'linux':
    pass
else:
    warnings.warn('The platform %s is not supported. '
                  'This may or may not work...' % sys.platform)

faf_mod = Extension('%s._fastqandfurious' % PACKAGENAME,
                    sources=['src/_fastqandfurious.c', ],
                    #depends=['src/.h'],
                    #include_dirs=['src',],
                    language='c',
                    extra_compile_args=(extra_compile_args +
                                        ['-O3', '-std=c99']))

setup(
    package_dir = {PACKAGENAME: 'src'}, 
    ext_modules = [faf_mod, ]
)
