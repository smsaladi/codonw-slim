import os
import glob

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize

import numpy as np

ext_files = glob.glob("codonw/codonwlib/src/*.c")
ext_files.extend(glob.glob("codonw/codonwlib/*.pyx"))

codonwlib = Extension(
    "codonw.codonwlib",
    ext_files,
    include_dirs=["codonw/codonwlib/include/", np.get_include()],
)

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='codonw-slim',
    version='1.5.0',
    license='GPLv2',
    author='Shyam Saladi',
    author_email='saladi@caltech.edu',
    url='https://github.com/smsaladi/codonw-slim',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=[
        'codonw',
    ],
    setup_requires=[
        'cython',
    ],
    install_requires=[
        'numpy',
        'pandas'
    ],
    tests_require = [
        'pytest',
        'biopython',
    ],
    test_suite="pytest",
    ext_modules=cythonize([codonwlib], language_level="3")
)
