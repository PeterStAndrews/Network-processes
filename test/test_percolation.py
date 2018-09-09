# test percolation for `Network-processes`
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

from network_processes import *
import unittest
import epyc

class PercolationTest(unittest.TestCase):
	'''Tests for `PERCOLATION` class in `percolation.py` using an
	epyc lab simulation environment. '''

	def setUp( self ):
		'''Set up the parameters.'''
		# set lab and test parameters
		self._lab = epyc.Lab()

		self._lab[PERCOLATION.T] = [0.6, 0.9]
		self._lab[PERCOLATION.N] = 5000
		self._lab[PERCOLATION.AVERAGE_K] = 5

		# repetitions at each point in the parameter space
		self._repetitions = 1

	def testEpidemic( self ):
		''''Test an epidemic occurs.'''
		# instance class
		e = PERCOLATION()
		# perform the experiments
		self._lab.runExperiment(epyc.RepeatedExperiment(e, self._repetitions))
		# extract the results
		rc = (self._lab.results())[0]
		# perform tests
		self.assertTrue(rc[epyc.Experiment.RESULTS]['occupied_fraction'] > 0)
		


