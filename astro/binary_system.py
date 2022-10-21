import numpy as np

GRAVITATION = 6.67428e-11


class BinarySystem:
    """Integrator for solving the two-body problem on Keplerian orbits.
    """
    def __init__(self, mass_1, radius_1, mass_2, radius_2, semi_major_axis, eccentricity, longitude_of_ascending_node,
                 inclination, argument_of_periapsis):
        """Initializes a gravitationally bound system with two masses.

        :param mass_1:                      Mass of the primary object in kilogram.
        :param radius_1:                    Radius of the primary object in meters.
        :param mass_2:                      Mass of the secondary object in kilogram.
        :param radius_2:                    Radius of the secondary object in meters.
        :param semi_major_axis:             Sum of the semi-distances between periapsis and apoapsis of both objects
                                            in meters.
        :param eccentricity:                Shape of the ellipses describing the orbits of both objects.
        :param longitude_of_ascending_node: Orientation of the semi-major axes within the reference plane in radians.
        :param inclination:                 Vertical tilt of the semi-major axes with respect to the reference plane
                                            in radians.
        :param argument_of_periapsis:       Orientation of the semi-major axes within the orbital plane in radians.
        """
        self.m_1 = mass_1
        self.m_2 = mass_2
        self.r_1 = radius_1
        self.r_2 = radius_2

        self.a = semi_major_axis
        self.e = eccentricity

        phi = longitude_of_ascending_node
        theta = inclination - 0.5 * np.pi
        psi = argument_of_periapsis

        self.p = np.sqrt(4 * np.pi**2 * self.a**3 / (GRAVITATION * (self.m_1 + self.m_2)))

        self.a_1 = self.a * self.m_2 / (self.m_1 + self.m_2)
        self.a_2 = self.a - self.a_1

        self.x_1 = self.a_1 * (1 + self.e) * (np.cos(phi) * np.cos(psi) - np.sin(phi) * np.cos(theta) * np.sin(psi))
        self.y_1 = self.a_1 * (1 + self.e) * (np.sin(phi) * np.cos(psi) + np.cos(phi) * np.cos(theta) * np.sin(psi))
        self.z_1 = self.a_1 * (1 + self.e) * np.sin(theta) * np.sin(psi)

        self.x_2 = -self.a_2 * (1 + self.e) * (np.cos(phi) * np.cos(psi) - np.sin(phi) * np.cos(theta) * np.sin(psi))
        self.y_2 = -self.a_2 * (1 + self.e) * (np.sin(phi) * np.cos(psi) + np.cos(phi) * np.cos(theta) * np.sin(psi))
        self.z_2 = -self.a_2 * (1 + self.e) * np.sin(theta) * np.sin(psi)

        v = np.sqrt(GRAVITATION * (self.m_1 + self.m_2) / self.a * (1 - self.e) / (1 + self.e))

        self.vx_1 = v * self.m_2 / (self.m_1 + self.m_2) * (np.cos(phi) * np.sin(psi)
                                                            + np.sin(phi) * np.cos(theta) * np.cos(psi))
        self.vy_1 = v * self.m_2 / (self.m_1 + self.m_2) * (np.sin(phi) * np.sin(psi)
                                                            - np.cos(phi) * np.cos(theta) * np.cos(psi))
        self.vz_1 = -v * self.m_2 / (self.m_1 + self.m_2) * np.sin(theta) * np.cos(psi)

        self.vx_2 = -v * self.m_1 / (self.m_1 + self.m_2) * (np.cos(phi) * np.sin(psi)
                                                             + np.sin(phi) * np.cos(theta) * np.cos(psi))
        self.vy_2 = -v * self.m_1 / (self.m_1 + self.m_2) * (np.sin(phi) * np.sin(psi)
                                                             - np.cos(phi) * np.cos(theta) * np.cos(psi))
        self.vz_2 = v * self.m_1 / (self.m_1 + self.m_2) * np.sin(theta) * np.cos(psi)

        self.orbit_1 = None
        self.orbit_2 = None

    def calc_stability_limit(self):
        """Calculates minimum distance from barycenter at which stable orbits are possible.

        return: Radius of minimum stable orbit in meters.
        """
        mu = self.m_2 / (self.m_1 + self.m_2)
        return self.a * (1.6 + 4.12 * mu + 5.1 * self.e - 4.27 * mu * self.e - 5.09 * mu**2 - 2.22 * self.e**2
                         + 4.61 * mu**2 * self.e**2)

    def calc_hill_sphere(self):
        """Calculates radius of Hill sphere where the gravitational forces of both objects are equal.

        return: Radius of Hill sphere in meters.
        """
        return self.a * (1 - self.e) * np.cbrt(self.m_2 / (3 * self.m_1))

    def calc_roche_limit(self):
        """Calculates Roche limit where tidal forces and gravitational self-attraction are equal.

        return: Radius of Roche limit in meters.
        """
        return self.r_2 * np.cbrt(2 * self.m_1 / self.m_2)

    def integrate_orbits(self, dt, t_max):
        """Integrates orbits using the Runge-Kutta method.

        :param dt:    Step size in seconds.
        :param t_max: Maximum integration time in seconds.
        """
        def func(y):
            r = np.sqrt((y[6] - y[0])**2 + (y[7] - y[1])**2 + (y[8] - y[2])**2)
            f = np.empty(12)
            f[0] = y[3]
            f[1] = y[4]
            f[2] = y[5]
            f[3] = -GRAVITATION * self.m_2 * (y[0] - y[6]) / r**3
            f[4] = -GRAVITATION * self.m_2 * (y[1] - y[7]) / r**3
            f[5] = -GRAVITATION * self.m_2 * (y[2] - y[8]) / r**3
            f[6] = y[9]
            f[7] = y[10]
            f[8] = y[11]
            f[9] = -GRAVITATION * self.m_1 * (y[6] - y[0]) / r**3
            f[10] = -GRAVITATION * self.m_1 * (y[7] - y[1]) / r**3
            f[11] = -GRAVITATION * self.m_1 * (y[8] - y[2]) / r**3
            return f

        n = np.ceil(t_max / dt).astype(np.uint32)
        t = 0.0

        self.orbit_1 = np.empty((n + 1, 7))
        self.orbit_2 = np.empty((n + 1, 7))

        y = np.array([self.x_1, self.y_1, self.z_1, self.vx_1, self.vy_1, self.vz_1,
                      self.x_2, self.y_2, self.z_2, self.vx_2, self.vy_2, self.vz_2])

        self.orbit_1[0] = [t, y[0], y[1], y[2], y[3], y[4], y[5]]
        self.orbit_2[0] = [t, y[6], y[7], y[8], y[9], y[10], y[11]]

        for i in range(1, n + 1):
            k_1 = dt * func(y)
            k_2 = dt * func(y + 0.5 * k_1)
            k_3 = dt * func(y + 0.5 * k_2)
            k_4 = dt * func(y + k_3)

            y += (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
            t += dt

            self.orbit_1[i] = [t, y[0], y[1], y[2], y[3], y[4], y[5]]
            self.orbit_2[i] = [t, y[6], y[7], y[8], y[9], y[10], y[11]]
