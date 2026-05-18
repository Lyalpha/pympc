import datetime
import json
import os
import tempfile
import unittest

import numpy as np
import numpy.testing as npt
from astropy.table import Table
from astropy.time import Time

import pympc
from pympc.pympc import (
    generate_xephem_catalogue,
    get_catalogue_status,
    XEPHEM_FILENAME_TEMPLATE,
)
from pympc.utils import get_observatory_data

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
TEST_MPCORB_XEPHEM = os.path.join(FILE_DIR, "resources", "mpcorb_xephem_head.csv")
TEST_MPCORB_XEPHEM_EXPECTED = os.path.join(
    FILE_DIR,
    "resources",
    "mpcorb_xephem_expected.csv",
)
TEST_MPCORB_XEPHEM_JSON = os.path.join(
    FILE_DIR,
    "resources",
    "mpcorb_xephem_test.json",
)


def _format_astorb_line(
    number="",
    name="",
    orbit_computer="E. Bowell",
    H="3.34",
    G=0.15,
    bv="",
    diameter="",
    taxon="",
    codes=(0, 0, 0, 0, 0, 0),
    arc=56959,
    num_obs=4750,
    epoch=(1996, 4, 27),
    mean_anomaly=80.477333,
    peri=71.802404,
    node=80.659857,
    inclination=10.600303,
    eccentricity=0.07604100,
    semimajor_axis=2.76788714,
    orbit_comp_date=(1996, 4, 14),
    ceu=0.0,
    ceu_rate=0.0,
    ceu_date=(1996, 4, 19),
    next_peu=(0.0, (1996, 5, 30)),
    max_peu=(0.0, (2004, 1, 11)),
    post_obs_peu=(0.0, (2004, 1, 11)),
):
    return (
        f"{str(number):>6} {str(name):<18} {str(orbit_computer):<15} {str(H):>5} {G:5.2f} "
        f"{str(bv):<4} {str(diameter):>5} {str(taxon):<4} "
        f"{codes[0]:4d}{codes[1]:4d}{codes[2]:4d}{codes[3]:4d}{codes[4]:4d}{codes[5]:4d} "
        f"{arc:5d}{num_obs:5d} {epoch[0]:04d}{epoch[1]:02d}{epoch[2]:02d} "
        f"{mean_anomaly:10.6f} {peri:10.6f} "
        f"{node:10.6f}{inclination:10.6f} {eccentricity:10.8f}{semimajor_axis:13.8f} "
        f"{orbit_comp_date[0]:04d}{orbit_comp_date[1]:02d}{orbit_comp_date[2]:02d} "
        f"{ceu:7.2f} {ceu_rate:8.2f} {ceu_date[0]:04d}{ceu_date[1]:02d}{ceu_date[2]:02d}"
        f" {next_peu[0]:7.2f} {next_peu[1][0]:04d}{next_peu[1][1]:02d}{next_peu[1][2]:02d}"
        f" {max_peu[0]:7.2f} {max_peu[1][0]:04d}{max_peu[1][1]:02d}{max_peu[1][2]:02d}"
        f" {post_obs_peu[0]:7.2f} {post_obs_peu[1][0]:04d}{post_obs_peu[1][1]:02d}{post_obs_peu[1][2]:02d}"
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
            if np.issubdtype(a.dtype, np.floating) or np.issubdtype(
                b.dtype, np.floating
            ):
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
                get_observatory_data(bad_obs)  # type: ignore

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
        )
        self._assert_tables_equal(ceres_result, self.ceres_result_geo)

    def test_iter_xephem_entries(self):
        with open(TEST_MPCORB_XEPHEM_EXPECTED, "r") as f:
            expected_lines = [line.strip() for line in f]

        iter_lines = list(pympc.pympc._iter_xephem_entries(TEST_MPCORB_XEPHEM_EXPECTED))
        self.assertEqual(iter_lines, expected_lines)

    def test_iter_xephem_chunks(self):
        with open(TEST_MPCORB_XEPHEM_EXPECTED, "r") as f:
            expected_lines = [line.strip() for line in f]

        for chunk_size, expected_chunk_lengths in ((3, [3, 3, 3, 1]), (4, [4, 4, 2])):
            with self.subTest(chunk_size=chunk_size):
                chunks = list(
                    pympc.pympc._iter_xephem_chunks(
                        TEST_MPCORB_XEPHEM_EXPECTED, chunk_size
                    )
                )
                self.assertEqual(
                    [len(chunk) for chunk in chunks], expected_chunk_lengths
                )
                self.assertTrue(all(chunk for chunk in chunks))
                self.assertEqual(
                    [line for chunk in chunks for line in chunk], expected_lines
                )

    def test_minor_planet_check_searchrad(self):
        # No result should be returned for a search radius of 0
        ceres_result = pympc.minor_planet_check(
            self.ceres_ra,
            self.ceres_dec,
            self.ceres_mjd,
            self.ceres_search_radius * 0,
            TEST_MPCORB_XEPHEM,
        )
        self.assertEqual(len(ceres_result), 0)

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
        )
        self.assertEqual(len(ceres_result), 0)

    def test_minor_planet_check_skips_entries_with_runtimeerror_on_mag(self):
        """Test that minor_planet_check gracefully handles nearly-parabolic bodies far from the Sun.

        This uses the same test body pyephem uses in their own tests.
        See: https://github.com/brandon-rhodes/pyephem/blob/master/pyephem/tests/test_ephem.py
        """
        parabolic_xephem_line = (
            "C/1980 Y1 (Bradfield),e,138.5850,115.3515,358.2941,945.0557,"
            "0.0000339,0.999725,359.9999,12/27.0/1980,2000,g  9.0,4.0"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(parabolic_xephem_line + "\n")
            temp_file = f.name

        try:
            with self.assertWarnsRegex(
                RuntimeWarning,
                "Could not compute properties for body",
            ):
                result = pympc.minor_planet_check(
                    self.ceres_ra,
                    self.ceres_dec,
                    self.ceres_mjd,
                    self.ceres_search_radius,
                    xephem_filepath=temp_file,
                    include_major_bodies=False,
                    chunk_size=0,
                )
            self.assertEqual(len(result), 0)
        finally:
            os.remove(temp_file)

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
        self.csv_file = None
        self.json_file = TEST_MPCORB_XEPHEM_JSON
        self.expected_csv_file = TEST_MPCORB_XEPHEM_EXPECTED

    def test_generate_mpcorb_xephem(self):
        self.csv_file = generate_xephem_catalogue(
            base_filepath=self.json_file,
            source="mpcorb",
        )

        self.assertTrue(os.path.exists(self.csv_file))

        with open(self.csv_file, "r") as f:
            generated_csv = f.read()
        with open(self.expected_csv_file, "r") as f:
            expected_csv = f.read()
        self.assertEqual(generated_csv, expected_csv)

    def tearDown(self):
        if self.csv_file is not None and os.path.exists(self.csv_file):
            os.remove(self.csv_file)


class TestGenerateASTORB(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.astorb_file = os.path.join(self.tmpdir.name, "astorb.dat")
        self.nea_file = os.path.join(self.tmpdir.name, "nea.json")
        self.csv_file = None

        astorb_lines = [
            _format_astorb_line(
                number="1",
                name="Ceres",
                H="3.34",
                G=0.12,
                bv="0.72",
                diameter="913.0",
                taxon="G?",
                arc=56959,
                num_obs=4750,
                epoch=(1996, 4, 27),
                mean_anomaly=80.477333,
                peri=71.802404,
                node=80.659857,
                inclination=10.600303,
                eccentricity=0.07604100,
                semimajor_axis=2.76788714,
                orbit_comp_date=(1996, 4, 14),
                ceu=0.03,
                ceu_rate=0.0,
                ceu_date=(1996, 4, 19),
                next_peu=(0.03, (1996, 5, 30)),
                max_peu=(0.03, (2004, 1, 11)),
                post_obs_peu=(0.03, (2004, 1, 11)),
            ),
            _format_astorb_line(
                number="",
                name="2024 AB",
                H="21.40",
                G=0.15,
                arc=120,
                num_obs=50,
                epoch=(2024, 9, 1),
                mean_anomaly=22.500000,
                peri=120.000000,
                node=88.000000,
                inclination=7.500000,
                eccentricity=0.12345678,
                semimajor_axis=1.98765432,
                orbit_comp_date=(2024, 8, 28),
                ceu=1.0,
                ceu_rate=0.5,
                ceu_date=(2024, 9, 1),
                next_peu=(2.0, (2024, 10, 1)),
                max_peu=(3.0, (2025, 1, 1)),
                post_obs_peu=(1.5, (2025, 1, 2)),
            ),
        ]
        with open(self.astorb_file, "w") as f:
            f.write("\n".join(astorb_lines) + "\n")

        nea_rows = [
            {
                "Number": "",
                "Name": "2024 AB",
                "Principal_desig": "2024 AB",
                "Epoch": Time("2024-09-02", scale="tt").jd,
                "M": 44.0,
                "Peri": 121.0,
                "Node": 89.0,
                "i": 8.0,
                "e": 0.2,
                "a": 1.8,
                "n": 360.0 / (365.25689 * 1.8**1.5),
                "H": 20.0,
                "G": 0.25,
            }
        ]
        with open(self.nea_file, "w") as f:
            json.dump(nea_rows, f)

    def test_generate_xephem_catalogue_astorb(self):
        self.csv_file = generate_xephem_catalogue(
            base_filepath=self.astorb_file,
            source="astorb",
        )
        self.assertTrue(os.path.exists(self.csv_file))
        with open(self.csv_file, "r") as f:
            generated_lines = [line.strip() for line in f if line.strip()]

        expected_epoch = Time("1996-04-27", scale="tt").utc.decimalyear
        expected_n = 360.0 / (365.25689 * 2.76788714**1.5)
        self.assertEqual(len(generated_lines), 2)
        self.assertEqual(
            generated_lines[0],
            f"Ceres,e,10.60030300,80.65985700,71.80240400,2.76788714,{expected_n:.8f},0.07604100,80.47733300,{expected_epoch:.8f},2000,3.34000000,0.12000000",
        )

    def test_generate_xephem_catalogue_astorb_with_nea_overlay(self):
        self.csv_file = generate_xephem_catalogue(
            base_filepath=self.astorb_file,
            nea_filepath=self.nea_file,
            source="astorb",
        )
        with open(self.csv_file, "r") as f:
            generated_lines = [line.strip() for line in f if line.strip()]

        self.assertEqual(len(generated_lines), 2)
        self.assertTrue(any(line.startswith("Ceres,e,") for line in generated_lines))
        self.assertTrue(
            any(
                line.startswith(
                    "2024 AB,e,8.00000000,89.00000000,121.00000000,1.80000000,"
                )
                for line in generated_lines
            )
        )
        self.assertFalse(
            any("2.00000000,0.12345678,22.50000000" in line for line in generated_lines)
        )

    def tearDown(self):
        if self.csv_file is not None and os.path.exists(self.csv_file):
            os.remove(self.csv_file)
        self.tmpdir.cleanup()


class TestCatalogueStatus(unittest.TestCase):
    def test_get_catalogue_status(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base_dir = os.path.join(tmpdir, "base")
            overlay_dir = os.path.join(tmpdir, "overlay")
            xephem_dir = os.path.join(tmpdir, "xephem")
            os.makedirs(base_dir, exist_ok=True)
            os.makedirs(overlay_dir, exist_ok=True)
            os.makedirs(xephem_dir, exist_ok=True)

            source = "astorb"
            base_file = os.path.join(base_dir, "astorb.dat.zst")
            nea_file = os.path.join(overlay_dir, "nea.json.zst")
            comets_file = os.path.join(overlay_dir, "cometels.json.zst")
            xephem_file = os.path.join(
                xephem_dir, XEPHEM_FILENAME_TEMPLATE.format(source=source)
            )

            for path in (base_file, nea_file, comets_file, xephem_file):
                with open(path, "w") as f:
                    f.write("test")

            status = get_catalogue_status(cat_dir=tmpdir, source=source)

            self.assertEqual(status["cache_dir"], tmpdir)
            self.assertEqual(status["source"], source)
            self.assertTrue(status["sources"]["base"]["exists"])
            self.assertTrue(status["sources"]["nea"]["exists"])
            self.assertTrue(status["sources"]["comets"]["exists"])
            self.assertTrue(status["xephem"]["exists"])
            self.assertEqual(status["xephem"]["path"], xephem_file)
