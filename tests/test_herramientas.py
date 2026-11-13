import unittest

import numpy as np
import sympy as sp

from sistemas_dinamicos import convolucion, escalon, integrar, routh


class HerramientasTest(unittest.TestCase):
    def test_escalon(self):
        np.testing.assert_array_equal(escalon([-1, 0, 1]), [0, 1, 1])

    def test_convolucion_causal(self):
        t = np.linspace(0, 2, 2001)
        y = convolucion(t, escalon(t), np.exp(-t))
        np.testing.assert_allclose(y, 1 - np.exp(-t), atol=0.0011)

    def test_convolucion_rechaza_malla_irregular(self):
        with self.assertRaises(ValueError):
            convolucion([0, 1, 3], [1, 1, 1], [1, 1, 1])

    def test_metodos_y_pasos_variables(self):
        t = np.array([0, 0.1, 0.3, 0.5, 1])
        f = lambda t, x: -x
        exacta = np.exp(-1)
        euler = integrar(f, t, [1], metodo="euler")[-1, 0]
        rk4 = integrar(f, t, [1], metodo="rk4")[-1, 0]
        self.assertLess(abs(rk4 - exacta), abs(euler - exacta))
        self.assertLess(abs(rk4 - exacta), 0.0002)

    def test_validacion_integracion(self):
        for tiempo, inicial, metodo in [
            ([0, 0], [1], "rk4"),
            ([0, 1], [], "rk4"),
            ([0, 1], [1], "otro"),
        ]:
            with self.subTest(tiempo=tiempo, inicial=inicial, metodo=metodo):
                with self.assertRaises(ValueError):
                    integrar(lambda t, x: -x, tiempo, inicial, metodo=metodo)
        with self.assertRaises(ValueError):
            integrar(lambda t, x: [1, 2], [0, 1], [1])

    def test_routh_estable_e_inestable(self):
        self.assertEqual(list(routh([1, 6, 11, 6])[:, 0]), [1, 6, 10, 6])
        self.assertEqual(
            list(routh([1, 2, -2, -1])[:, 0]), [1, 2, -sp.Rational(3, 2), -1]
        )
        with self.assertRaises(ValueError):
            routh([1, 0, 1])
        with self.assertRaises(ValueError):
            routh([1, sp.oo, 2])


if __name__ == "__main__":
    unittest.main()
