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
            config = yaml.safe_load(file)

        self.mass_stars_in_sun = config[star_system].get('mass stars', [1.0, 1.0])
        self.radius_stars_in_sun = config[star_system].get('radius stars', [1.0, 1.0])
        self.luminosity_stars_in_sun = config[star_system].get('luminosity stars', [1.0, 1.0])
        self.temperature_stars_in_kelvin = config[star_system].get('temperature stars', [5780, 5780])
        self.semi_major_axis_binary_star_in_au = config[star_system].get('semi-major axis stars', 1.0)
        self.eccentricity_binary_star = config[star_system].get('eccentricity stars', 0.0)
        self.inclination_binary_star_in_deg = config[star_system].get('inclination stars', 0.0)
        self.rotation_binary_star_in_deg = config[star_system].get('rotation stars', 0.0)

        self.mass_planets_in_jupiter = config[star_system].get('mass planets', [])
        self.radius_planets_in_jupiter = config[star_system].get('radius planets', [])
        self.semi_major_axis_planets_in_au = config[star_system].get('semi-major axis planets', [])
        self.eccentricity_planets = config[star_system].get('eccentricity planets', [])
        self.inclination_planets_in_deg = config[star_system].get('inclination planets', [])
        self.rotation_planets_in_deg = config[star_system].get('rotation planets', [])

        self.mass_satellites_in_earth = config[star_system].get('mass satellites', [[]])
        self.radius_satellites_in_earth = config[star_system].get('radius satellites', [[]])
        self.semi_major_axis_satellites_in_au = config[star_system].get('semi-major axis satellites', [[]])
        self.eccentricity_satellites = config[star_system].get('eccentricity satellites', [[]])
        self.inclination_satellites_in_deg = config[star_system].get('inclination satellites', [[]])
        self.rotation_satellites_in_deg = config[star_system].get('rotation satellites', [[]])
