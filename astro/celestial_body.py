import numpy as np

GRAVITATION = 6.67428e-11


class CelestialBody:
    """Spherical object that follows Kepler's laws.
    """
    def __init__(self, mass, radius):
        """Sets properties of celestial body.

        :param mass:   Mass of the object.
        :param radius: Mean radius of the object.
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

        :param mass:        Mass of the object.
        :param radius:      Mean radius of the object.
        :param luminosity:  Bolometric luminosity of the object.
        :param temperature: Effective temperature of the object.
        """
        super().__init__(mass, radius)

        self.luminosity = luminosity
        self.temperature = temperature
