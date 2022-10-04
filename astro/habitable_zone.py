import numpy as np

AU_TO_METER = 149597870700
SUN_TEMPERATURE = 5780
SUN_TO_WATT = 3.828e26


class HabitableZone:
    """Circumbinary habitable zone.
    """
    def __init__(self, star_1, star_2):
        """Initializes luminosities for calculating habitability limits.

        :param star_1: Primary star.
        :param star_2: Secondary star.
        """
        t_1 = star_1.temperature
        t_2 = star_2.temperature
        l_1 = star_1.luminosity / SUN_TO_WATT
        l_2 = star_2.luminosity / SUN_TO_WATT

        self.l_1 = l_1 / self.calc_effective_stellar_flux(t_1)
        self.l_2 = l_2 / self.calc_effective_stellar_flux(t_2)

    @staticmethod
    def calc_effective_stellar_flux(t_eff):
        """Calculates insolation thresholds at different atmospheric collapse limits.

        :param t_eff: Effective temperature of the star in Kelvin.

        :return Effective stellar flux at recent Venus, runaway greenhouse, maximum greenhouse and early Mars limit.
        """
        t = t_eff - SUN_TEMPERATURE
        s = np.empty(4)

        s[0] = 1.776 + 2.136e-4 * t + 2.533e-8 * t**2 - 1.332e-11 * t**3 - 3.097e-15 * t**4
        s[1] = 1.107 + 1.332e-4 * t + 1.580e-8 * t**2 - 8.308e-12 * t**3 - 1.931e-15 * t**4
        s[2] = 0.356 + 6.171e-5 * t + 1.698e-9 * t**2 - 3.198e-12 * t**3 - 5.575e-16 * t**4
        s[3] = 0.320 + 5.547e-5 * t + 1.526e-9 * t**2 - 2.874e-12 * t**3 - 5.011e-16 * t**4

        return s

    def calc_radiative_habitable_limits(self, p_1, p_2):
        """Calculates boundaries of the radiative habitable zone for different habitability models.

        :param p_1: Cartesian coordinate of primary star at apastron.
        :param p_2: Cartesian coordinate of secondary star at apastron.

        :return Radii of inner and outer radiative habitable limits in meters.
        """
        d = 0.5 * np.linalg.norm(np.subtract(p_1, p_2)) / AU_TO_METER
        r_in, r_out, _, _ = self.solve_differential_equation(d, 0)
        return np.abs(np.hstack((r_in[:2], r_out[2:]))) * AU_TO_METER

    def solve_differential_equation(self, a, phi):
        """Solves quartic equation for calculating habitability limits.

        For detailed information, see https://iopscience.iop.org/article/10.1088/0004-637X/780/1/14/pdf.

        :param a:   Semi-distance between both stars in astronomical units.
        :param phi: True anomaly of hypothetical orbiting object on habitability limit contour in radians.

        return: Four possible solutions of equation ordered by descending relevance.
        """
        a_0 = a**4 - a**2 * (self.l_1 + self.l_2)
        a_1 = 2 * a * np.cos(phi) * (self.l_1 - self.l_2)
        a_2 = 2 * a**2 * (1 - 2 * np.cos(phi)**2) - self.l_1 - self.l_2

        b_0 = 4 * a_0 * a_2 - a_1**2
        b_1 = -4 * a_0
        b_2 = -a_2

        q = b_1 / 3 - b_2**2 / 9
        r = -b_0 / 2 + b_1 * b_2 / 6 - b_2**3 / 27
        s = np.cbrt(r + np.sqrt(q**3 + r**2))
        t = np.cbrt(r - np.sqrt(q**3 + r**2))

        y_1 = -b_2 / 3 + (s + t)

        c = np.sqrt(-2 * a**2 * (1 - 2 * np.cos(phi)**2) + self.l_1 + self.l_2 + y_1)
        d = np.sqrt(self.l_1 + self.l_2 + 4 * a * (self.l_1 - self.l_2) / c * np.cos(phi)
                    - 2 * a**2 * (1 - 2 * np.cos(phi)**2) - y_1)
        e = np.sqrt(self.l_1 + self.l_2 - 4 * a * (self.l_1 - self.l_2) / c * np.cos(phi)
                    - 2 * a ** 2 * (1 - 2 * np.cos(phi) ** 2) - y_1)

        return -c/2 - d/2, -c/2 + d/2, c/2 - e/2, c/2 + e/2
