import csv
import os
import unittest

import pympc

TEST_MPCORB_XEPHEM = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "mpcorb_xephem_head.csv"
)


class TestPyMPC(unittest.TestCase):
    def setUp(self):
        self.ceres_ra = 61.78375
        self.ceres_dec = 21.945
        self.ceres_search_radius = 5
        self.ceres_mjd = 59640.0
        self.ceres_result = pympc.pympc._to_astropy_table(
            [
                [
                    "Ceres",
                    61.78350721167082,
                    21.945118347316505,
                    0.9158510091022737,
                    8.82,
                    "Ceres,e,10.58769000,80.26858000,73.63703000,2.76604310,0.21424745,0.07850100,291.37563000,2022.05479452,2000,3.56000000,0.15000000",
                ]
            ]
        )

    # def test_update_catalogue(self):
    #     catalogue = pympc.update_catalogue(cat_dir=".")
    #     with open(catalogue, "r") as f:
    #         csv.Sniffer().sniff(f.read(), delimiters=",")
    #     os.remove(catalogue)

    def test_minor_planet_check(self):
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            chunk_size=0,
        )
        self.assertEqual(ceres_result, self.ceres_result)

    def test_minor_planet_check_multiprocess(self):
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius,
            TEST_MPCORB_XEPHEM,
            chunk_size=2e4,
        )
        self.assertEqual(ceres_result, self.ceres_result)

    def test_minor_planet_check_searchrad(self):
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
