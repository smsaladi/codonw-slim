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

setup(
    name='CodonW',
    version='1.0',
    author='Shyam Saladi',
    author_email='saladi@caltech.edu',
    url='https://github.com/smsaladi/codonw-slim',
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
