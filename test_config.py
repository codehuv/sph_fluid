import unittest
from config import Config, N, SIM_W, BOTTOM, DAM, DAM_BREAK, G, SPACING, K, K_NEAR, REST_DENSITY, R, SIGMA, MAX_VEL, WALL_DAMP, VEL_DAMP, GRID_CELL_SIZE

class TestConfig(unittest.TestCase):
    def test_return_config(self):
        config_instance = Config()
        params = config_instance.return_config()

        self.assertEqual(len(params), 16, "Config.return_config() should return 16 parameters.")

        self.assertEqual(params[0], N)
        self.assertEqual(params[1], SIM_W)
        self.assertEqual(params[2], BOTTOM)
        self.assertEqual(params[3], DAM)
        self.assertEqual(params[4], DAM_BREAK)
        self.assertAlmostEqual(params[5], G)
        self.assertAlmostEqual(params[6], SPACING)
        self.assertAlmostEqual(params[7], K)
        self.assertAlmostEqual(params[8], K_NEAR)
        self.assertAlmostEqual(params[9], REST_DENSITY)
        self.assertAlmostEqual(params[10], R)
        self.assertAlmostEqual(params[11], SIGMA)
        self.assertAlmostEqual(params[12], MAX_VEL)
        self.assertAlmostEqual(params[13], WALL_DAMP)
        self.assertAlmostEqual(params[14], VEL_DAMP)
        self.assertAlmostEqual(params[15], GRID_CELL_SIZE)

        # Check some types
        self.assertIsInstance(N, int)
        self.assertIsInstance(SIM_W, int) # or float, depending on use
        self.assertIsInstance(G, float)
        self.assertIsInstance(R, float)

if __name__ == '__main__':
    unittest.main()
