import unittest
import numpy as np
from zibell.simulate import data, simul_rep
from zibell.estimation import estime_rep
from zibell.bell import bell_pmf, bell_rvs

class TestSimulation(unittest.TestCase):
    
    def test_data_output_shapes(self):
        n = 100
        b = np.array([-0.5, 1.2])
        g = -1.1
        result = data(n, b, g)
        self.assertEqual(result['Y'].shape[0], n)
        self.assertEqual(result['X'].shape, (2, n))
        self.assertEqual(result['Z'].shape[0], n)

    def test_simul_rep_dimensions(self):
        n = 100
        nbrep = 10
        b = np.array([-0.5, 1.2])
        g = 0.2
        result = simul_rep(n, nbrep, b, g)
        self.assertEqual(result['My'].shape, (nbrep, n))
        self.assertEqual(result['Mx'].shape, (2, n, nbrep))

    def test_estimation_runs(self):
        n = 50
        b = np.array([-0.5, 1.2])
        g = -1.1
        d = data(n, b, g)
        est = estime_rep(d['Y'], d['X'], d['Z'], alpha=0.1)
        self.assertEqual(len(est), 3)
        self.assertTrue(np.all(np.isfinite(est)))  # Pas de NaN ou inf

    def test_bell_functions(self):
        mu = 3
        k = 2
        prob = bell_pmf(k, mu)
        self.assertTrue(0 <= prob <= 1)
        sample = bell_rvs(mu, size=10)
        self.assertEqual(len(sample), 10)

if __name__ == '__main__':
    unittest.main()
