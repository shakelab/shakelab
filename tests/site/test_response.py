# =============================================================================
#
# Copyright (C) 2010-2017 GEM Foundation
#
# This file is part of the OpenQuake's Site Response Toolkit (OQ-SRTK)
#
# OQ-SRTK is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# OQ-SRTK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# Author: Valerio Poggi
#
# =============================================================================

import unittest
import numpy as np
import numpy.testing as npt

from shakelab.site.response import impedance_amplification
from shakelab.site.response import sh_transfer_function


# =============================================================================

class ImpedanceAmplificationTestCase(unittest.TestCase):
    """
    Testing the Impedance Contrast (IC) method to
    compute amplification
    """

    def check_amplification(self,
                            top_vs,
                            top_dn,
                            ref_vs,
                            ref_dn,
                            inc_angle,
                            expected_result,
                            tolerance=0.):

        computed_result = impedance_amplification(top_vs,
                                                  top_dn,
                                                  ref_vs,
                                                  ref_dn,
                                                  inc_angle)

        if not isinstance(computed_result, np.ndarray):
            self.assertAlmostEquals(computed_result,
                                    expected_result,
                                    delta=tolerance)
        else:
            for cr, er, in zip(computed_result, expected_result):
                self.assertAlmostEquals(cr,
                                        er,
                                        delta=tolerance)

    def test_no_interface(self):
        """
        Case of a model with no seismic contrast
        """

        self.check_amplification(200., 1900., 200., 1900., 10., 1.)

    def test_one_interface_vertical(self):
        """
        Case of a model with a single seismic impedance contrast
        """

        self.check_amplification(200., 1900., 1500., 2500., 0., 3.14, 0.01)

    def test_one_interface_angle(self):
        """
        Case of a model with a single seismic impedance contrast
        and non-vertical incidence angle
        """

        self.check_amplification(200., 1900., 1500., 2500., 45., 2.65, 0.01)

    def test_multiple_interfaces_vertical(self):
        """
        Case of a model with an array of seismic impedance contrasts
        """

        self.check_amplification(np.array([200., 800., 2000.]),
                                 np.array([1900., 2100., 2500.]),
                                 2000.,
                                 2500.,
                                 0.,
                                 np.array([3.63, 1.72, 1.]),
                                 0.01)


# =============================================================================

class ShTansferFunctionTestCase(unittest.TestCase):
    """
    Test for the calculation of the SH-wave transfer
    function at arbitrary depth

    Results are verified against those obtained from the PSVQ code
    of F.J. Sanchez Sesma, modified by Roberto Paolucci
    """

    def check_amplification(self,
                            test_file,
                            hl,
                            vs,
                            dn,
                            qs,
                            inc_ang,
                            depth,
                            decimal=2.):

        freq, ampf = self.read_psvq_file(test_file)
        freq = np.array(freq)

        disp = sh_transfer_function(freq, hl, vs, dn, qs, inc_ang, depth)

        npt.assert_almost_equal(ampf,
                                np.abs(disp[0]/2),
                                decimal=decimal)

    def read_psvq_file(self, input_file):

        freq = []
        ampf = []
        with open(input_file, 'r') as f:
            for line in f:
                value = line.strip().split()
                freq.append(float(value[0]))
                ampf.append(float(value[1]))
            f.close()

        return freq, ampf

    def test_2_layer_model_free_surface(self):
        """
        Case with two layers
        """

        self.check_amplification('tests/site/data/psvq/2_layer_model/amp.i0.txt',
                                 np.array([100., 0.]),
                                 np.array([200., 1200.]),
                                 np.array([1900., 2500.]),
                                 np.array([10., 100.]),
                                 0.,
                                 0.)

    def test_2_layer_model_1st_interface(self):
        """
        Case with two layers and calculation at the 1st interface
        """

        self.check_amplification('tests/site/data/psvq/2_layer_model/amp.i1.txt',
                                 np.array([100., 0.]),
                                 np.array([200., 1200.]),
                                 np.array([1900., 2500.]),
                                 np.array([10., 100.]),
                                 0.,
                                 100.)

    def test_3_layer_model_free_surface(self):
        """
        Case with three layers
        """

        self.check_amplification('tests/site/data/psvq/3_layer_model/amp.i0.txt',
                                 np.array([10., 50., 0]),
                                 np.array([200., 500., 1200.]),
                                 np.array([1900., 2100., 2500.]),
                                 np.array([10., 20., 100.]),
                                 0.,
                                 0.)

    def test_3_layer_model_1st_interface(self):
        """
        Case with three layers and calculation at the 1st interface
        """

        self.check_amplification('tests/site/data/psvq/3_layer_model/amp.i1.txt',
                                 np.array([10., 50., 0]),
                                 np.array([200., 500., 1200.]),
                                 np.array([1900., 2100., 2500.]),
                                 np.array([10., 20., 100.]),
                                 0.,
                                 10.)

    def test_3_layer_model_2nd_interface(self):
        """
        Case with three layers and calculation at the 2st interface
        """

        self.check_amplification('tests/site/data/psvq/3_layer_model/amp.i2.txt',
                                 np.array([10., 50., 0]),
                                 np.array([200., 500., 1200.]),
                                 np.array([1900., 2100., 2500.]),
                                 np.array([10., 20., 100.]),
                                 0.,
                                 60.)

if __name__ == '__main__':
    unittest.main()
