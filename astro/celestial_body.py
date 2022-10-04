import numpy as np

GRAVITATION = 6.67428e-11


class CelestialBody:
    """Spherical object that follows Kepler's laws.
    """
    def __init__(self, mass, radius):
        """Sets properties of celestial body.

        :param mass:   Mass in kilogram.
        :param radius: Mean radius in meters.
        """
        self.mass = mass
        self.radius = radius

        self.density = 3 * self.mass / (4 * np.pi * self.radius**3)
        self.gravity = GRAVITATION * self.mass / self.radius**2
        self.escape_velocity = np.sqrt(2 * GRAVITATION * self.mass / self.radius)


class Star(CelestialBody):
    """Celestial body that has sufficient mass to trigger nuclear fusion.
    """
    def __init__(self, mass, radius, luminosity, temperature):
        """Sets properties of star.

        :param mass:        Mass in kilogram.
        :param radius:      Mean radius in meters.
        :param luminosity:  Bolometric luminosity in Watt.
        :param temperature: Effective temperature in Kelvin.
        """
        super().__init__(mass, radius)

        self.luminosity = luminosity
        self.temperature = temperature
