#!/usr/bin/env python3
"""Input-only spectrum contract tests; run on Linux, without a solver/build."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET


HERE = Path(__file__).resolve().parent
PREPARE = HERE / "prepare.py"
TABLE = HERE.parent / "cbc1971_theta42_spectrum_3d.dat"
BASE = ["--length", ".6", "--width", ".6", "--h", ".3", "--modes", "4",
        "--hom-modes", "32", "--dt", ".0001", "--duration", ".002",
        "--record-start", ".001", "--cutoff", "100", "--stations", "0,.3",
        "--buffer", ".15", "--probe-grid", "2,2"]


class PrepareSpectrum(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="hit-prepare-spectrum-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def create(self, name, options=(), expect_success=True):
        directory = self.root / name
        result = subprocess.run([sys.executable, str(PREPARE), "--output", str(directory)]
                                + BASE + list(options), capture_output=True, text=True)
        if not expect_success:
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((directory / "pilot.xml").exists(), result.stderr)
            return result.stderr
        self.assertEqual(result.returncode, 0, result.stderr)
        return (directory, ET.parse(directory / "pilot.xml").getroot(),
                json.loads((directory / "manifest.json").read_text()))

    def test_cbc_default_matches_explicit_selection(self):
        first, root, manifest = self.create("default")
        second, _, _ = self.create("explicit", ["--spectrum", "cbc"])
        self.assertEqual((first / "pilot.xml").read_bytes(), (second / "pilot.xml").read_bytes())
        self.assertEqual((first / TABLE.name).read_bytes(), TABLE.read_bytes())
        self.assertEqual(root.find("./HITINFLOW/REGION/SPECTRUM").attrib,
                         {"TYPE": "Tabulated", "FILE": TABLE.name, "INTERPOLATION": "LogLinear",
                          "LOW": "K4", "HIGH": "Gaussian", "NORMALIZATION": "Preserve"})
        self.assertEqual(root.find("./HITINFLOW/REGION/CUTOFF").attrib,
                         {"KMAX": "100", "TAPERSTART": "80"})
        self.assertEqual(manifest["physics"], {"mean_velocity_m_s": 12.7, "kinvis_m2_s": 1.5e-5,
            "spectrum": "CBC1971 Table 3 theta=42 shell-3D; upstream only", "normalization": "Preserve",
            "interpolation": "LogLinear", "low_extension": "K4", "high_extension": "Gaussian"})
        self.assertEqual(manifest["sha256"][TABLE.name], hashlib.sha256(TABLE.read_bytes()).hexdigest())

    def test_von_karman_preserves_amplitude_without_unused_table(self):
        directory, root, manifest = self.create("vk", ["--spectrum", "von-karman", "--vk-length",
            ".031", "--vk-amplitude", "2.5", "--taper", "80"])
        spectrum = root.find("./HITINFLOW/REGION/SPECTRUM")
        self.assertEqual(set(spectrum.attrib), {"TYPE", "LENGTH", "AMPLITUDE", "NORMALIZATION"})
        self.assertEqual((spectrum.get("TYPE"), spectrum.get("NORMALIZATION")), ("VonKarman", "Preserve"))
        self.assertEqual(float(spectrum.get("LENGTH")), .031)
        self.assertEqual(float(spectrum.get("AMPLITUDE")), 2.5)
        self.assertEqual(root.find("./HITINFLOW/REGION/CUTOFF").attrib,
                         {"KMAX": "100", "TAPERSTART": "80"})
        self.assertFalse((directory / TABLE.name).exists())
        self.assertNotIn(TABLE.name, manifest["sha256"])
        physics = manifest["physics"]
        self.assertEqual((physics["spectrum_type"], physics["normalization"]), ("VonKarman", "Preserve"))
        self.assertEqual((physics["vk_length_m"], physics["vk_amplitude"]), (.031, 2.5))
        self.assertIn("distinct", physics["length_note"])
        self.assertIn("divergent", physics["uncut_dissipation_note"])
        self.assertIn("no retained-band renormalization", physics["retained_band_note"])
        self.assertIn("performs no fit", physics["parameter_provenance"]["source"])

    def test_zero_amplitude_is_valid(self):
        _, root, manifest = self.create("zero", ["--spectrum", "von-karman", "--vk-length",
            ".031", "--vk-amplitude", "0", "--taper", "80"])
        self.assertEqual(root.find("./HITINFLOW/REGION/SPECTRUM").get("AMPLITUDE"), "0")
        self.assertEqual(manifest["physics"]["vk_amplitude"], 0)

    def test_matrix_dealiasing_fails_before_creating_an_invalid_session(self):
        error = self.create("invalid-matrix", ["--backend", "Matrix", "--dealias"], expect_success=False)
        self.assertIn("homogeneous dealiasing requires --backend FFTW", error)
        self.assertFalse((self.root / "invalid-matrix").exists())
        _, root, _ = self.create("matrix", ["--backend", "Matrix"])
        self.assertIsNone(root.find("./CONDITIONS/SOLVERINFO/I[@PROPERTY='USEFFT']"))
        self.assertIsNone(root.find("./CONDITIONS/SOLVERINFO/I[@PROPERTY='DEALIASING']"))

    def test_invalid_or_unused_spectrum_options_fail_before_writing(self):
        vk = ["--spectrum", "von-karman", "--taper", "80"]
        cases = [["--vk-length", ".03"], ["--vk-amplitude", "1"],
                 vk, vk + ["--vk-length", ".03"], vk + ["--vk-amplitude", "1"],
                 ["--spectrum", "von-karman", "--vk-length", ".03", "--vk-amplitude", "1"]]
        cases += [vk + ["--vk-length", value, "--vk-amplitude", "1"]
                  for value in ("0", "-1", "nan", "inf")]
        cases += [vk + ["--vk-length", ".03", "--vk-amplitude", value]
                  for value in ("-1", "nan", "inf")]
        for index, options in enumerate(cases):
            with self.subTest(options=options):
                self.create(f"invalid-{index}", options, expect_success=False)


if __name__ == "__main__":
    unittest.main()
