import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation

from astro.binary_system import BinarySystem
from utils.parameter_config import ParameterConfig


AU_TO_METER = 149597870700
DAY_TO_SEC = 86400
DEG_TO_RAD = np.pi / 180
EARTH_TO_KG = 5.97237e24
EARTH_TO_METER = 6371000
JUPITER_TO_KG = 1.8982e27
JUPITER_TO_METER = 69911000
SUN_TO_KG = 1.98847e30


# Script parameters
star_system = "Kepler-16"
planet_index = 0
resolution_in_deg = 2


def simulate(num, planet, satellites, hill_sphere=None):
    """Simulates movement of satellites on their orbits around the host planet.

    :param num:         Ordinal number of planetary system.
    :param planet:      Host planet.
    :param satellites:  Satellites with orbit histories to be simulated.
    :param hill_sphere: Radius of planet's Hill sphere.
    """
    fig, ax = plt.subplots(num=f"Planet {num+1}")
    fig.set_facecolor('black')

    ax.set_aspect('equal', 'datalim')
    ax.set_axis_off()

    if hill_sphere is not None:
        ax.add_patch(Circle((0, 0), hill_sphere, fill=False, edgecolor='red', linestyle='--'))

    ax.add_patch(Circle((0, 0), planet.r_2 * 5, color='blue'))

    position_satellites = []
    for satellite in satellites:
        period = np.flatnonzero(satellite.orbit_2[:, 0] >= satellite.p)[0] + 1
        ax.plot(satellite.orbit_2[:period, 1], satellite.orbit_2[:period, 2], color='white')
        init_position = Circle((satellite.orbit_2[0, 1], satellite.orbit_2[0, 2]), satellite.r_2 * 50, color='gray')
        position_satellites.append(ax.add_patch(init_position))

    time = ax.text(0.0, 0.97, "Orbital period: 0 days", color='white', transform=ax.transAxes)

    def update(i):
        for n, satellite in enumerate(satellites):
            position_satellites[n].set_center((satellite.orbit_2[i, 1], satellite.orbit_2[i, 2]))
        time.set_text(f"Orbital period: {satellites[0].orbit_2[i, 0] / DAY_TO_SEC:.2f} days")
        return *position_satellites, time

    ani = FuncAnimation(fig, update, frames=len(satellites[0].orbit_2), interval=30, blit=True)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    config = ParameterConfig(star_system)

    planet = BinarySystem(
        sum(config.mass_stars_in_sun) * SUN_TO_KG, None,
        config.mass_planets_in_jupiter[planet_index] * JUPITER_TO_KG,
        config.radius_planets_in_jupiter[planet_index] * JUPITER_TO_METER,
        config.semi_major_axis_planets_in_au[planet_index] * AU_TO_METER,
        config.eccentricity_planets[planet_index],
        config.inclination_planets_in_deg[planet_index] * DEG_TO_RAD,
        config.rotation_planets_in_deg[planet_index] * DEG_TO_RAD
    )

    hill_sphere = planet.calc_hill_sphere()

    num_satellites = len(config.mass_satellites_in_earth[planet_index])
    satellites = []

    # planet has to contain at least one satellite
    assert num_satellites > 0

    for n in range(num_satellites):
        satellite = BinarySystem(
            planet.m_2, planet.r_2,
            config.mass_satellites_in_earth[planet_index][n] * EARTH_TO_KG,
            config.radius_satellites_in_earth[planet_index][n] * EARTH_TO_METER,
            config.semi_major_axis_satellites_in_au[planet_index][n] * AU_TO_METER,
            config.eccentricity_satellites[planet_index][n],
            config.inclination_satellites_in_deg[planet_index][n] * DEG_TO_RAD,
            config.rotation_satellites_in_deg[planet_index][n] * DEG_TO_RAD
        )

        # satellite has to be located beyond Roche limit but within its planet's Hill sphere
        assert hill_sphere > satellite.a_2 > satellite.calc_roche_limit()

        satellites.append(satellite)

    step_size = min([satellite.p for satellite in satellites]) * resolution_in_deg / 360
    simulation_time = max([satellite.p for satellite in satellites])

    [satellite.integrate_orbits(step_size, simulation_time) for satellite in satellites]

    simulate(planet_index, planet, satellites, hill_sphere)
