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

        self.semi_major_axis_binary_star_in_au = config.get('semi-major axis', 0.2)
        self.eccentricity_binary_star = config.get('eccentricity', 0.0)
        self.longitude_of_ascending_node_binary_star_in_deg = config.get('longitude of ascending node', 0.0)
        self.inclination_binary_star_in_deg = config.get('inclination', 90.0)
        self.argument_of_periapsis_binary_star_in_deg = config.get('argument of periapsis', 0.0)

        stars = config.get('stars', [])
        self.mass_stars_in_sun = self.aggregate(stars, 'mass', 1.0)
        self.radius_stars_in_sun = self.aggregate(stars, 'radius', 1.0)
        self.luminosity_stars_in_sun = self.aggregate(stars, 'luminosity', 1.0)
        self.temperature_stars_in_kelvin = self.aggregate(stars, 'temperature', 5780)

        planets = config.get('planets', [])
        self.mass_planets_in_jupiter = self.aggregate(planets, 'mass', 1.0)
        self.radius_planets_in_jupiter = self.aggregate(planets, 'radius', 1.0)
        self.semi_major_axis_planets_in_au = self.aggregate(planets, 'semi-major axis', 1.0)
        self.eccentricity_planets = self.aggregate(planets, 'eccentricity', 0.0)
        self.longitude_of_ascending_node_planets_in_deg = self.aggregate(planets, 'longitude of ascending node', 0.0)
        self.inclination_planets_in_deg = self.aggregate(planets, 'inclination', 90.0)
        self.argument_of_periapsis_planets_in_deg = self.aggregate(planets, 'argument of periapsis', 0.0)

        satellites = [planet.get('satellites', []) for planet in planets]
        self.mass_satellites_in_earth = self.aggregate(satellites, 'mass', 0.1)
        self.radius_satellites_in_earth = self.aggregate(satellites, 'radius', 0.5)
        self.semi_major_axis_satellites_in_au = self.aggregate(satellites, 'semi-major axis', 0.01)
        self.eccentricity_satellites = self.aggregate(satellites, 'eccentricity', 0.0)
        self.longitude_of_ascending_node_satellites_in_deg = self.aggregate(
            satellites, 'longitude of ascending node', 0.0)
        self.inclination_satellites_in_deg = self.aggregate(satellites, 'inclination', 90.0)
        self.argument_of_periapsis_satellites_in_deg = self.aggregate(satellites, 'argument of periapsis', 0.0)

    @staticmethod
    def aggregate(config, key, default):
        """Searches configuration recursively and aggregates all values that match key.

        :param config:  Configuration to be searched.
        :param key:     Key to be matched.
        :param default: Default value if key is not available.
        """
        aggregations = []

        for element in config:
            if isinstance(element, list):
                aggregations.append(ParameterConfig.aggregate(element, key, default))
            else:
                aggregations.append(element.get(key, default))

        return aggregations
