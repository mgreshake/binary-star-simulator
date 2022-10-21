import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from matplotlib.animation import FuncAnimation

from astro.binary_system import BinarySystem
from astro.celestial_body import Star
from astro.habitable_zone import HabitableZone
from utils.parameter_config import ParameterConfig


AU_TO_METER = 149597870700
DAY_TO_SEC = 86400
DEG_TO_RAD = np.pi / 180
JUPITER_TO_KG = 1.8982e27
JUPITER_TO_METER = 69911000
SUN_TO_KG = 1.98847e30
SUN_TO_METER = 695700000
SUN_TO_WATT = 3.828e26


# Script parameters
star_system = "Kepler-16"
is_habitability_simulated = True
are_planets_simulated = True
resolution_in_deg = 2


def simulate(star_system, binary_star, stability_limit=None, habitability_limit=None, planets=()):
    """Simulates movement of celestial objects on their orbits within a binary star system.

    :param star_system:        Name of binary star system.
    :param binary_star:        Binary star with orbit history to be simulated.
    :param stability_limit:    Radius of minimum stable orbit around binary star.
    :param habitability_limit: Radii of different inner and outer habitable limits.
    :param planets:            Planets with orbit histories to be simulated.
    """
    fig, ax = plt.subplots(num=star_system)
    fig.set_facecolor('black')

    ax.set_aspect('equal', 'datalim')
    ax.set_axis_off()

    if habitability_limit is not None:
        ax.add_patch(Wedge((0, 0), habitability_limit[3], 0, 360, width=habitability_limit[3] - habitability_limit[0],
                           facecolor='green', edgecolor=None, alpha=0.5))
        ax.add_patch(Wedge((0, 0), habitability_limit[2], 0, 360, width=habitability_limit[2] - habitability_limit[1],
                           facecolor='green', edgecolor=None, alpha=0.5))

    if stability_limit is not None:
        ax.add_patch(Circle((0, 0), stability_limit, fill=False, edgecolor='red', linestyle='--'))

    period = np.flatnonzero(binary_star.orbit_1[:, 0] >= binary_star.p - 1e-07)[0] + 1
    ax.plot(binary_star.orbit_1[:period, 1], binary_star.orbit_1[:period, 2], color='white')
    ax.plot(binary_star.orbit_2[:period, 1], binary_star.orbit_2[:period, 2], color='white')

    radius_star_1 = binary_star.r_1 * 5
    radius_star_2 = binary_star.r_2 * 5

    init_position_star_1 = Circle((binary_star.orbit_1[0, 1], binary_star.orbit_1[0, 2]), radius_star_1, color='yellow')
    init_position_star_2 = Circle((binary_star.orbit_2[0, 1], binary_star.orbit_2[0, 2]), radius_star_2, color='red')

    position_star_1 = ax.add_patch(init_position_star_1)
    position_star_2 = ax.add_patch(init_position_star_2)

    position_planets = []
    for planet in planets:
        period = np.flatnonzero(planet.orbit_2[:, 0] >= planet.p - 1e-07)[0] + 1
        ax.plot(planet.orbit_2[:period, 1], planet.orbit_2[:period, 2], color='white')
        init_position = Circle((planet.orbit_2[0, 1], planet.orbit_2[0, 2]), planet.r_2 * 30, color='blue')
        position_planets.append(ax.add_patch(init_position))

    time = ax.text(0.0, 0.97, "Orbital period: 0 days", color='white', transform=ax.transAxes)

    def update(i):
        position_star_1.set_center((binary_star.orbit_1[i, 1], binary_star.orbit_1[i, 2]))
        position_star_2.set_center((binary_star.orbit_2[i, 1], binary_star.orbit_2[i, 2]))
        for n, planet in enumerate(planets):
            position_planets[n].set_center((planet.orbit_2[i, 1], planet.orbit_2[i, 2]))
        time.set_text(f"Orbital period: {binary_star.orbit_1[i, 0] / DAY_TO_SEC:.2f} days")
        return position_star_1, position_star_2, *position_planets, time

    ani = FuncAnimation(fig, update, frames=len(binary_star.orbit_1), interval=30, blit=True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    config = ParameterConfig(star_system)

    star_1 = Star(
        config.mass_stars_in_sun[0] * SUN_TO_KG,
        config.radius_stars_in_sun[0] * SUN_TO_METER,
        config.luminosity_stars_in_sun[0] * SUN_TO_WATT,
        config.temperature_stars_in_kelvin[0]
    )

    star_2 = Star(
        config.mass_stars_in_sun[1] * SUN_TO_KG,
        config.radius_stars_in_sun[1] * SUN_TO_METER,
        config.luminosity_stars_in_sun[1] * SUN_TO_WATT,
        config.temperature_stars_in_kelvin[1]
    )

    binary_star = BinarySystem(
        star_1.mass, star_1.radius,
        star_2.mass, star_2.radius,
        config.semi_major_axis_binary_star_in_au * AU_TO_METER,
        config.eccentricity_binary_star,
        config.longitude_of_ascending_node_binary_star_in_deg * DEG_TO_RAD,
        config.inclination_binary_star_in_deg * DEG_TO_RAD,
        config.argument_of_periapsis_binary_star_in_deg * DEG_TO_RAD
    )

    stability_limit = binary_star.calc_stability_limit()

    if is_habitability_simulated:
        habitable_zone = HabitableZone(star_1, star_2)

        habitability_limit = habitable_zone.calc_radiative_habitable_limits(
            (binary_star.x_1, binary_star.y_1),
            (binary_star.x_2, binary_star.y_2)
        )

        # inner habitability limits have to be located completely inside outer habitability limits
        assert (habitability_limit[:-1] <= habitability_limit[1:]).all()
    else:
        habitability_limit = None

    if are_planets_simulated:
        # orbital period of binary star has to be at least 7 days
        assert binary_star.p > 604800

        num_planets = len(config.semi_major_axis_planets_in_au)
        planets = []

        for n in range(num_planets):
            planet = BinarySystem(
                star_1.mass + star_2.mass, None,
                config.mass_planets_in_jupiter[n] * JUPITER_TO_KG,
                config.radius_planets_in_jupiter[n] * JUPITER_TO_METER,
                config.semi_major_axis_planets_in_au[n] * AU_TO_METER,
                config.eccentricity_planets[n],
                config.longitude_of_ascending_node_planets_in_deg[n] * DEG_TO_RAD,
                config.inclination_planets_in_deg[n] * DEG_TO_RAD,
                config.argument_of_periapsis_planets_in_deg[n] * DEG_TO_RAD
            )

            # distance from planet to barycenter has to be greater than minimum stable orbit around binary star
            assert planet.a_2 > stability_limit

            planets.append(planet)
    else:
        planets = []

    step_size = binary_star.p * resolution_in_deg / 360
    simulation_time = max([planet.p for planet in planets]) if len(planets) > 0 else binary_star.p

    binary_star.integrate_orbits(step_size, simulation_time)
    [planet.integrate_orbits(step_size, simulation_time) for planet in planets]

    simulate(star_system, binary_star, stability_limit, habitability_limit, planets)
