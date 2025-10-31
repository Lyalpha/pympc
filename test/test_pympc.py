import datetime
import os
import unittest

import numpy as np
import numpy.testing as npt
from astropy.table import Table

import pympc
from pympc.pympc import generate_mpcorb_xephem
from pympc.utils import get_observatory_data

TEST_MPCORB_XEPHEM = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "resources", "mpcorb_xephem_head.csv"
)


class TestPyMPC(unittest.TestCase):

    def _assert_tables_equal(self, result, expected, rtol=1e-6, atol=0):
        """
        Compare two astropy.table.Table objects in tests.

        - Verifies column names and row count.
        - For float columns uses assert_allclose (or decimal-based comparison).
        - For masked columns ensures masks match and compares filled values.
        - For non-float columns uses exact array equality.
        """

        if not (isinstance(result, Table) and isinstance(expected, Table)):
            self.fail("Both result and expected must be astropy.table.Table instances")

        # schema
        self.assertListEqual(list(result.colnames), list(expected.colnames))
        self.assertEqual(len(result), len(expected))

        for name in result.colnames:
            a = result[name].data
            b = expected[name].data

            # normalize masked arrays
            if np.ma.isMaskedArray(a) or np.ma.isMaskedArray(b):
                a = np.ma.asanyarray(a)
                b = np.ma.asanyarray(b)
                # masks must match
                if not np.array_equal(a.mask, b.mask):
                    self.fail(f"Mask differs for column {name!s}")
                a = a.filled(np.nan)
                b = b.filled(np.nan)

            # floats: tolerant compare
            if np.issubdtype(a.dtype, np.floating) or np.issubdtype(b.dtype, np.floating):
                npt.assert_allclose(a, b, rtol=rtol, atol=atol)
            else:
                # exact compare for ints/strings/objects
                if not np.array_equal(a, b):
                    self.fail(f"Column {name!s} differs")

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

        self.moon_ra = 350.725037
        self.moon_dec = -6.13578
        self.moon_mjd = 60481.147089
        self.moon_search_radius = 5
        self.moon_result_topo = pympc.pympc._to_astropy_table(
            [
                [
                    "Hyperion",
                    350.7246696793632,
                    -6.135791897494748,
                    1.315472890149291,
                    14.1,
                    "Hyperion,P",
                ],
            ]
        )

    def test_observatory_data_retrieval(self):
        # Check the observatory data is as expected
        geocentric_tuple = (0.0, 0.0, 0.0)
        for obs in (500, "500", "Geocentric", geocentric_tuple):
            self.assertEqual(get_observatory_data(obs), geocentric_tuple)
            
        lapalma_tuple = (342.1176, 0.87764, 0.47847)
        for obs in (950, "950", "La Palma", lapalma_tuple):
            self.assertEqual(get_observatory_data(obs), lapalma_tuple)

        # Check a non-existent observatory returns ValueError
        for bad_obs in (1234.5, "ZZZ", "NonExistent Observatory", (1, 2), (1, 2, 3, 4)):
            with self.assertRaises(ValueError):
                get_observatory_data(bad_obs)

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
        self._assert_tables_equal(ceres_result, self.ceres_result_geo)

    def test_minor_planet_check_multiprocess(self):
        # The same result should be returned for a multiprocessing call
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            chunk_size=20000,
        )
        self._assert_tables_equal(ceres_result, self.ceres_result_geo)

    def test_minor_planet_check_searchrad(self):
        # No result should be returned for a search radius of 0
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius * 0,
            TEST_MPCORB_XEPHEM,
            chunk_size=20000,
        )
        self.assertEqual(ceres_result, [])

        # Check a 100% sky coverage search radius returns all minor objects
        all_results = pympc.minor_planet_check(
            0,
            0,
            self.ceres_mjd,
            3600 * 180,
            TEST_MPCORB_XEPHEM,
            include_major_bodies=False,
            chunk_size=0,
        )
        self.assertEqual(len(all_results), 10)

        # Check a 100% sky coverage search radius returns all major objects
        all_results = pympc.minor_planet_check(
            0,
            0,
            self.ceres_mjd,
            3600 * 180,
            TEST_MPCORB_XEPHEM,
            include_minor_bodies=False,
            chunk_size=0,
        )
        self.assertEqual(len(all_results), 27)

    def test_minor_planet_check_maxmag(self):
        # No results should be returned for a maxmag of 0
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            max_mag=0,
            chunk_size=20000,
        )
        self.assertEqual(ceres_result, [])

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
            self._assert_tables_equal(ceres_result, self.ceres_result_geo)
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
            self._assert_tables_equal(ceres_result, self.ceres_result_topo)

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
                rhocosphi=0,
                rhosinphi=0,
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
                rhocosphi=0,
                rhosinphi=0,
            )
            self.assertEqual((ra_topo, dec_topo), (ra_geo, dec_geo))
        # A distance of 0 should return a ZeroDivisionError
        dist_au = 0
        with self.assertRaises(ZeroDivisionError):
            _ = pympc.pympc.equitorial_geocentric_to_topocentric(
                ra_geo=0,
                dec_geo=0,
                dist_au=dist_au,
                epoch=epoch,
                longitude=0,
                rhocosphi=0,
                rhosinphi=0,
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
                rhocosphi=0.5,
                rhosinphi=0.5,
            )
            npt.assert_allclose(_topo_coo_ret, _coo_topo, rtol=1e-6, atol=0)

    def test_major_moon_check(self):
        moon_result = pympc.minor_planet_check(
            self.moon_ra,
            self.moon_dec,
            self.moon_mjd,
            self.moon_search_radius,
            include_minor_bodies=False,
            observatory=950,
            chunk_size=0,
        )
        self._assert_tables_equal(moon_result, self.moon_result_topo)


class TestGenerateMPCORB(unittest.TestCase):
    def setUp(self):
        script_dir = os.path.dirname(os.path.abspath(__file__))

        self.csv_file = None
        self.json_file = os.path.join(
            script_dir, "resources", "mpcorb_xephem_test.json"
        )
        self.expected_csv_file = os.path.join(
            script_dir, "resources", "mpcorb_xephem_expected.csv"
        )

    def test_generate_mpcorb_xephem(self):
        self.csv_file = generate_mpcorb_xephem(self.json_file)

        self.assertTrue(os.path.exists(self.csv_file))

        with open(self.csv_file, "r") as f:
            generated_csv = f.read()
        with open(self.expected_csv_file, "r") as f:
            expected_csv = f.read()
        self.assertEqual(generated_csv, expected_csv)

    def tearDown(self):
        if self.csv_file is not None and os.path.exists(self.csv_file):
            os.remove(self.csv_file)
