# Setup for `Network-processes`
#
# Copyright (C) 2017 Peter Mann
# 
# This file is part of `Network_processes`, for epidemic network 
# analytical results using Python.
#
# `Network_processes` is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# `Network_processes` is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with `Network_processes`. If not, see <http://www.gnu.org/licenses/gpl.html>.

from setuptools import setup

with open('README.rst') as f:
    longDescription = f.read()

setup(name = 'NetworkProcesses',
      version = '0.1.0',
      description = 'Network process simulation in Python',
      long_description = longDescription,
      url = 'https://github.com/PeterStAndrews/NetworkProcesses',
      author = 'Peter Mann',
      author_email = 'pm78@st-andrews.ac.uk',
      license = 'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      classifiers = [ 'Development Status :: Beta',
                      'Intended Audience :: Science/Research',
                      'Intended Audience :: Developers',
                      'Programming Language :: Python :: 2.7',
                      'Topic :: Scientific/Engineering' ],
      packages = [ 'NetworkProcesses' ],
      zip_safe = True,
      install_requires = [ "numpy", "networkx", "epyc", "scipy"])


