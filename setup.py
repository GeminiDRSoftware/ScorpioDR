#!/usr/bin/env python

"""
Setup script for scorpiodr

Usage: pip install [-e] .

"""

import os

from setuptools import setup, find_packages, Extension

# try:
#     from Cython.Build import cythonize
# except ImportError:
#     use_cython = False
# else:
#     use_cython = True

VERSION = '0.1.0'

PACKAGENAME = 'scorpiodr'
PACKAGES = find_packages()

# EXTENSIONS
# suffix = 'pyx' if use_cython else 'c'
# EXTENSIONS = [
#     Extension("gempy.library.cython_utils",
#               [os.path.join('gempy', 'library', 'cython_utils.' + suffix)])
# ]
# if use_cython:
#     EXTENSIONS = cythonize(EXTENSIONS)
EXTENSIONS = []

setup(name='scorpiodr',
      version=VERSION,
      description='Gemini Data Processing Python Package',
      author='Gemini Data Processing Software Group',
      author_email='sus_inquiries@gemini.edu',
      url='http://www.gemini.edu',
      maintainer='Science User Support Department',
      license='BSD',
      zip_safe=False,
      packages=PACKAGES,
      include_package_data=True,
      scripts=None,
      ext_modules=EXTENSIONS,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Gemini Ops',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: POSIX :: Linux',
          'Operating System :: Linux :: CentOS',
          'Operating System :: MacOS :: MacOS X',
          'Programming Language :: Python',
          'Topic :: Gemini',
          'Topic :: Data Reduction',
          'Topic :: Scientific/Engineering :: Astronomy',
      ],
      install_requires=[
          'dragons>=3.2',
      ],
      extras_require={
          'docs': ['sphinx', 'sphinx_rtd_theme'],
          'test': ['pytest', 'pytest_dragons>=1.0.0', 'coverage', 'objgraph'],
      },
      project_urls={
          'Issue Tracker': 'https://github.com/GeminiDRSoftware/ScorpioDR',
          'Documentation': 'https://dragons.readthedocs.io/',
      },
      # keywords=['astronomy', 'astrophysics', 'science', 'gemini'],
      python_requires='>=3.10',
      )
