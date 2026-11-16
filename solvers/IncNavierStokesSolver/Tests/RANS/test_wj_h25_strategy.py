# SPDX-License-Identifier: MIT
import math,unittest
from study_wj_h25 import measure_matches

class ProductionStrategy(unittest.TestCase):
    def test_derived_physical_measure_roundoff_is_not_state_payload_equality(self):
        self.assertTrue(measure_matches(25*(20-.0875-.15),494.0625))
        self.assertTrue(measure_matches(math.nextafter(497.8125,math.inf),497.8125))
        self.assertFalse(measure_matches(25*(20-.15),494.0625))
        self.assertFalse(measure_matches(494.0625+1e-7,494.0625))

    def test_invalid_or_wrong_measures_fail(self):
        for x in (math.nan,math.inf,-math.inf,3.75):self.assertFalse(measure_matches(x,497.8125))

if __name__=='__main__':unittest.main()
