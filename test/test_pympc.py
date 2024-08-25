import datetime
import os
import unittest

import numpy as np
from astropy.table import Table

import pympc
from pympc.utils import get_observatory_data

TEST_MPCORB_XEPHEM = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "mpcorb_xephem_head.csv"
)


class TestPyMPC(unittest.TestCase):
    def _roundedAssertEqual(self, result, expected, precision=8):
        """
        A custom comparison function for comparing results that removes astrophysically insignificant differences
        in the float values.
        """
        if isinstance(result, Table) and isinstance(expected, Table):
            result.round(precision)
            expected.round(precision)

        self.assertEqual(result, expected)

    def setUp(self):
        self.ceres_ra = 61.78375
        self.ceres_dec = 21.945
        self.ceres_mjd = 59640.0
        self.ceres_search_radius = 5
        self.ceres_result_geo = pympc.pympc._to_astropy_table(
            [
                [
                    "Ceres",
                    61.783502361039375,
                    21.945116761568123,
                    0.9276118611600402,
                    8.82,
                    "Ceres,e,10.58769000,80.26858000,73.63703000,2.76604310,0.21424745,0.07850100,291.37563000,2022.05479452,2000,3.56000000,0.15000000",
                ]
            ]
        )
        self.ceres_result_topo = pympc.pympc._to_astropy_table(
            [
                [
                    "Ceres",
                    61.78262858754154,
                    21.94475359752442,
                    3.848205908280207,
                    8.82,
                    "Ceres,e,10.58769000,80.26858000,73.63703000,2.76604310,0.21424745,0.07850100,291.37563000,2022.05479452,2000,3.56000000,0.15000000",
                ]
            ]
        )

    def test_observatory_data_retrieval(self):
        # Check the observatory data is as expected
        obs_data = get_observatory_data("500")
        self.assertEqual(obs_data, (0.0, 0.0, 0.0))
        obs_data = get_observatory_data("950")
        self.assertEqual(obs_data, (342.1176, 0.87764, 0.47847))

        # Check a non-existent observatory returns ValueError
        with self.assertRaises(ValueError):
            get_observatory_data("NONEXISTENT")

    def test_minor_planet_check(self):
        # Check the result is as expected
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            chunk_size=0,
        )
        self._roundedAssertEqual(ceres_result, self.ceres_result_geo)

    def test_minor_planet_check_multiprocess(self):
        # The same result should be returned for a multiprocessing call
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            chunk_size=2e4,
        )
        self._roundedAssertEqual(ceres_result, self.ceres_result_geo)

    def test_minor_planet_check_searchrad(self):
        # No result should be returned for a search radius of 0
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius * 0,
            TEST_MPCORB_XEPHEM,
            chunk_size=2e4,
        )
        self.assertEqual(ceres_result, None)

        # Check a 100% sky coverage search radius returns all objects
        all_results = pympc.minor_planet_check(
            0, 0, self.ceres_mjd, 3600 * 180, TEST_MPCORB_XEPHEM, chunk_size=0
        )
        self.assertEqual(len(all_results), 10)

    def test_minor_planet_check_maxmag(self):
        # No results should be returned for a maxmag of 0
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            max_mag=0,
            chunk_size=2e4,
        )
        self.assertEqual(ceres_result, None)

    def test_minor_planet_check_observatory(self):
        # Passing an observatory should return the correct results, regardless of
        # observatory format
        # Test geocentric coordinates are as expected
        for observatory in ("500", 500, (0, 0, 0)):
            ceres_result = pympc.minor_planet_check(
                self.ceres_ra,
                self.ceres_dec,
                self.ceres_mjd,
                self.ceres_search_radius,
                TEST_MPCORB_XEPHEM,
                observatory=observatory,
            )
            self.assertEqual(ceres_result, self.ceres_result_geo)
        # Test topocentric La Palma results are as expected
        for observatory in ("950", 950, (342.1176, 0.87764, +0.47847)):
            ceres_result = pympc.minor_planet_check(
                self.ceres_ra,
                self.ceres_dec,
                self.ceres_mjd,
                self.ceres_search_radius,
                TEST_MPCORB_XEPHEM,
                observatory=observatory,
            )
            self._roundedAssertEqual(ceres_result, self.ceres_result_topo)

        # Test a non-existent observatory returns ValueError
        with self.assertRaises(ValueError):
            pympc.minor_planet_check(
                self.ceres_ra,
                self.ceres_dec,
                self.ceres_mjd,
                self.ceres_search_radius,
                TEST_MPCORB_XEPHEM,
                observatory="NONEXISTENT",
            )

    def test_equatorial_geocentric_to_topocentric(self):
        coo_geo = [
            (0, 0),
            (0, np.pi / 2),
            (0, -np.pi / 2),
            (np.pi, 0),
            (np.pi, np.pi / 2),
            (np.pi, -np.pi / 2),
            (-np.pi, 0),
            (-np.pi, np.pi / 2),
            (-np.pi, -np.pi / 2),
            (2 * np.pi, 0),
            (2 * np.pi, np.pi / 2),
            (2 * np.pi, -np.pi / 2),
        ]
        epoch = datetime.datetime(2000, 1, 1)
        # An observer at the center of the earth should return the same coordinates for any position
        for ra_geo, dec_geo in coo_geo:
            ra_topo, dec_topo = pympc.pympc.equitorial_geocentric_to_topocentric(
                ra_geo=ra_geo,
                dec_geo=dec_geo,
                dist_au=1,
                epoch=epoch,
                longitude=0,
                rho_cos_phi=0,
                rho_sin_phi=0,
            )
            self.assertEqual((ra_topo, dec_topo), (ra_geo, dec_geo))
        # An observer at the center of the earth should return the same coordinates for any longitude
        ra_geo = dec_geo = 0
        for longitude in (0, 90, 180, 270, 360, -90, -180, -270):
            ra_topo, dec_topo = pympc.pympc.equitorial_geocentric_to_topocentric(
                ra_geo=0,
                dec_geo=0,
                dist_au=1,
                epoch=epoch,
                longitude=longitude,
                rho_cos_phi=0,
                rho_sin_phi=0,
            )
            self.assertEqual((ra_topo, dec_topo), (ra_geo, dec_geo))
        # A distance of 0 should return a ZeroDivisionError
        dist_au = 0
        with self.assertRaises(ZeroDivisionError):
            ra_topo, dec_topo = pympc.pympc.equitorial_geocentric_to_topocentric(
                ra_geo=0,
                dec_geo=0,
                dist_au=dist_au,
                epoch=epoch,
                longitude=0,
                rho_cos_phi=0,
                rho_sin_phi=0,
            )
        # An observer not at the center of the earth should return different coordinates for different longitudes
        coo_topo = [
            (-2.099587455215055e-05, -2.1317525541019556e-05),
            (-1.3968614631936522, 1.5707750087362562),
            (-1.3968614631936522, -1.5707750096451367),
            (3.141613649619262, -2.1317682831124705e-05),
            (4.886323843980277, -1.570817644853537),
            (4.886323843980277, 1.5708176439446566),
            (-3.1415716575603243, -2.1317682831124705e-05),
            (-1.3968614631993106, -1.570817644853537),
            (-1.3968614631993106, 1.5708176439446566),
            (6.283164311305034, -2.1317525541019556e-05),
            (4.8863238439859344, 1.5707750087362562),
            (4.8863238439859344, -1.5707750096451367),
        ]
        for _coo_geo, _coo_topo in zip(coo_geo, coo_topo):
            ra_topo, dec_topo = _coo_geo
            _topo_coo_ret = pympc.pympc.equitorial_geocentric_to_topocentric(
                ra_geo=ra_topo,
                dec_geo=dec_topo,
                dist_au=1,
                epoch=epoch,
                longitude=0,
                rho_cos_phi=0.5,
                rho_sin_phi=0.5,
            )
            self._roundedAssertEqual(_topo_coo_ret, _coo_topo)
