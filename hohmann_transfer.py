import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np

ORBITAL_CONSTANT = 398600
EARTH_RADIUS = 6378


def circular_orbit_velocity(radius_km):
    return math.sqrt(ORBITAL_CONSTANT / radius_km)


def elliptical_orbit_velocity(radius_apogee_km, radius_perigee_km):
    e2 = (radius_apogee_km - radius_perigee_km) / \
        (radius_apogee_km + radius_perigee_km)
    h2 = math.sqrt((radius_perigee_km / (1 / (1 + e2))) * ORBITAL_CONSTANT)
    v_perigee = h2 / radius_perigee_km
    v_apogee = h2 / radius_apogee_km
    return (v_perigee, v_apogee, e2)


def calculate_delta_v_manuevers(orbit_1_radius_km, orbit_2_radius_km):
    if (orbit_1_radius_km > orbit_2_radius_km):
        orbit_1_radius_km, orbit_2_radius_km = orbit_2_radius_km, orbit_1_radius_km
    v1 = circular_orbit_velocity(orbit_1_radius_km)
    v_perigee, v_apogee, e2 = elliptical_orbit_velocity(
        orbit_2_radius_km, orbit_1_radius_km)
    v3 = circular_orbit_velocity(orbit_2_radius_km)
    return (v1, v_perigee, v_apogee, v3, e2)


def calculate_total_required_delta_v(v1, v_perigee, v_apogee, v3):
    dv_1 = v_perigee - v1
    dv_2 = v3 - v_apogee
    dv_total = dv_1 + dv_2
    return dv_total


def plot_sphere():
    # TODO add ability for other orbit burn points.
    r = EARTH_RADIUS
    u, v = np.mgrid[0:2*np.pi:500j, 0:np.pi:500j]
    x = r * np.cos(u) * np.sin(v)
    y = r * np.sin(u) * np.sin(v)
    z = r * np.cos(v)
    return (x, y, z)


def plot_delta_v(orbit_1_radius_km, orbit_2_radius_km, projection_3d=False):
    figure, axes = plt.subplots(subplot_kw=dict(
        projection='3d' if projection_3d else None))

    max_radius = max(orbit_1_radius_km, orbit_2_radius_km)
    max_radius_scaled = max_radius + (max_radius * 0.2)
    axes.set_ylim(-max_radius_scaled, max_radius_scaled)
    axes.set_xlim(-max_radius_scaled, max_radius_scaled)

    if projection_3d:
        axes.set_zlim(-max_radius_scaled, max_radius_scaled)

    if (orbit_1_radius_km > orbit_2_radius_km):
        orbit_1_radius_km, orbit_2_radius_km = orbit_2_radius_km, orbit_1_radius_km

    circular_orbit_1 = plt.Circle((0, 0), orbit_1_radius_km, fill=False)
    circular_orbit_2 = plt.Circle(
        (0, 0), orbit_2_radius_km, fill=False)

    v1, v_perigee, v_apogee, v3, e2 = calculate_delta_v_manuevers(
        orbit_1_radius_km, orbit_2_radius_km)
    total_dv = calculate_total_required_delta_v(
        v1, v_perigee, v_apogee, v3)
    parameters = (f"Total dV = {total_dv}")
    print(parameters)

    center = ((-orbit_2_radius_km) + orbit_1_radius_km) / 2
    semi_minor_axis = (orbit_2_radius_km+orbit_1_radius_km) * \
        math.sqrt(1 - e2**2)

    ellipse = Ellipse((center, 0),
                      orbit_1_radius_km + orbit_2_radius_km, semi_minor_axis, color=(1, 0, 0), fill=False)

    for patch in [circular_orbit_1, circular_orbit_2, ellipse]:
        p = axes.add_patch(patch)
        if projection_3d:
            art3d.pathpatch_2d_to_3d(p, z=0, zdir="x")

    if projection_3d:
        x, y, z = plot_sphere()  # TODO add other orbit burn points.
        axes.plot(x, y, z, alpha=0.5)
    else:
        earth = plt.Circle((0, 0), EARTH_RADIUS)
        plt.plot(
            orbit_1_radius_km, 0, color=(0, 1, 0), marker="*")
        plt.plot(
            -orbit_2_radius_km, 0, color=(0, 1, 0), marker="*")
        axes.add_patch(earth)

    plt.show()


def main():
    plot_delta_v(orbit_1_radius_km=60000+EARTH_RADIUS,
                 orbit_2_radius_km=1000+EARTH_RADIUS, projection_3d=True)


if __name__ == "__main__":
    main()
