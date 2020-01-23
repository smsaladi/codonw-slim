import glob

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

codonwlib = Extension(
    "codonw.codonwlib",
    [*glob.glob("codonw/codonwlib/src/*.c"), "codonw/codonwlib/codonw.pyx"],
    include_dirs=["codonw/codonwlib/include/", np.get_include()],
)

setup(
    name='CodonW',
    version='1.0',
    author='Shyam Saladi',
    author_email='saladi@caltech.edu',
    url='https://github.com/smsaladi/codonw-barebones',
    packages=[
        'codonw',
    ],
    ext_modules=cythonize([codonwlib], language_level="3")
)
