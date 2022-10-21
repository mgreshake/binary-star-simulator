import yaml


class ParameterConfig:
    """Container that holds parameters for a binary star system.
    """
    def __init__(self, star_system, config_path="star_catalogue.yml"):
        """Reads parameters from configuration file.

        :param star_system: Binary star system to be parameterized.
        :param config_path: Path to configuration file.
        """
        with open(config_path, 'r') as file:
            content = yaml.safe_load(file)

        config = content[star_system]

        self.mass_stars_in_sun = config.get('mass stars', [1.0, 1.0])
        self.radius_stars_in_sun = config.get('radius stars', [1.0, 1.0])
        self.luminosity_stars_in_sun = config.get('luminosity stars', [1.0, 1.0])
        self.temperature_stars_in_kelvin = config.get('temperature stars', [5780, 5780])
        self.semi_major_axis_binary_star_in_au = config.get('semi-major axis stars', 0.2)
        self.eccentricity_binary_star = config.get('eccentricity stars', 0.0)
        self.longitude_of_ascending_node_binary_star_in_deg = config.get('longitude of ascending node stars', 0.0)
        self.inclination_binary_star_in_deg = config.get('inclination stars', 90.0)
        self.argument_of_periapsis_binary_star_in_deg = config.get('argument of periapsis stars', 0.0)

        self.semi_major_axis_planets_in_au = config.get('semi-major axis planets', [])
        planets = range(len(self.semi_major_axis_planets_in_au))
        self.mass_planets_in_jupiter = config.get('mass planets', [1.0 for _ in planets])
        self.radius_planets_in_jupiter = config.get('radius planets', [1.0 for _ in planets])
        self.eccentricity_planets = config.get('eccentricity planets', [0.0 for _ in planets])
        self.longitude_of_ascending_node_planets_in_deg = config.get(
            'longitude of ascending node planets', [0.0 for _ in planets])
        self.inclination_planets_in_deg = config.get('inclination planets', [90.0 for _ in planets])
        self.argument_of_periapsis_planets_in_deg = config.get('argument of periapsis planets', [0.0 for _ in planets])

        self.semi_major_axis_satellites_in_au = config.get('semi-major axis satellites', [[] for _ in planets])
        satellites = [range(len(self.semi_major_axis_satellites_in_au[p])) for p in planets]
        self.mass_satellites_in_earth = config.get(
            'mass satellites', [[0.1 for _ in satellites[p]] for p in planets])
        self.radius_satellites_in_earth = config.get(
            'radius satellites', [[0.5 for _ in satellites[p]] for p in planets])
        self.eccentricity_satellites = config.get(
            'eccentricity satellites', [[0.0 for _ in satellites[p]] for p in planets])
        self.longitude_of_ascending_node_satellites_in_deg = config.get(
            'longitude of ascending node satellites', [[0.0 for _ in satellites[p]] for p in planets])
        self.inclination_satellites_in_deg = config.get(
            'inclination satellites', [[90.0 for _ in satellites[p]] for p in planets])
        self.argument_of_periapsis_satellites_in_deg = config.get(
            'argument of periapsis satellites', [[0.0 for _ in satellites[p]] for p in planets])
