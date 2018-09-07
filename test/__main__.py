# run test suite for `Network-processes`
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

import unittest
from .test_network import *
from .test_percolation import *
from .test_gfs import *
from .test_sto import *
from .test_hmf import *
from .test_add_del import *


# initialise the tests
networkSuite = unittest.TestLoader().loadTestsFromTestCase(NetworkTest)
gfsSuite = unittest.TestLoader().loadTestsFromTestCase(GFsTest)
percolationSuite = unittest.TestLoader().loadTestsFromTestCase(PercolationTest)
stoSuite = unittest.TestLoader().loadTestsFromTestCase(STOTest)
hmfSuite = unittest.TestLoader().loadTestsFromTestCase(HMFTest)
addition_deletionSuite = unittest.TestLoader().loadTestsFromTestCase(addition_deletionTest)

# add tests to the test suite
suite = unittest.TestSuite([ networkSuite,
							 gfsSuite,
							 percolationSuite,
							 stoSuite,
							 hmfSuite,
							 addition_deletionSuite ] )

# run the tests
if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)












