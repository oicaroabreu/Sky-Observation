import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from skyfield.api import Star, load
from skyfield.api import N,S,E,W, wgs84
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection

import numpy as np

celestial_bodies = [
    'MERCURY BARYCENTER',
    'VENUS BARYCENTER',
    'MARS BARYCENTER',
    'JUPITER BARYCENTER',
    'SATURN BARYCENTER',
    'URANUS BARYCENTER',
    'NEPTUNE BARYCENTER',
    'PLUTO BARYCENTER',
    'SUN',
    'MERCURY',
    'VENUS',
    'MOON',
]


ts = load.timescale()
"""
Modify Date
Example of custom date: 

t = ts.utc(1993, 5, 15, 12, 15) 
"""
t = ts.now()
planets = load('de440s.bsp')
earth = planets['EARTH']

"""
Modify Location
Use latitude and longitude of the observation point on earth
"""
latitude = 16.7130
longitude = 49.2473
observer = wgs84.latlon(latitude * N, longitude * W)
observer_on_earth = earth + observer

ras = []
decs = []

"""
Celestial Bodies from Solar System
"""
for body in celestial_bodies:
    planet = planets[body]
    planet_pos = observer_on_earth.at(t).observe(planet)
    app = planet_pos.apparent()
    ra, dec, _ = app.radec()
    print(f"\n{body}:")
    print(ra.dstr(warn=False).split('deg')[0].strip())
    print(dec.dstr(warn=False).split('deg')[0].strip())
    ras.append(int(ra.dstr(warn=False).split('deg')[0].strip()))
    decs.append(int(dec.dstr(warn=False).split('deg')[0].strip()))

"""
Stars and Constelations
"""

with load.open(hipparcos.URL) as f:
    stars = hipparcos.load_dataframe(f)

url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
       '/skycultures/modern_st/constellationship.fab')

with load.open(url) as f:
    constellations = stellarium.parse_constellations(f)

edges = [edge for name, edges in constellations for edge in edges]
edges_star1 = [star1 for star1, star2 in edges]
edges_star2 = [star2 for star1, star2 in edges]

alt, az, _ = observer_on_earth.at(t).altaz()
center = observer.at(t).from_altaz(alt_degrees=alt.degrees, az_degrees=az.degrees)
projection = build_stereographic_projection(center)

star_positions = observer_on_earth.at(t).observe(Star.from_dataframe(stars))
stars['x'], stars['y'] = projection(star_positions)

print(stars)

limiting_magnitude = 4.40

bright_stars = (stars.magnitude <= limiting_magnitude)

magnitude = stars['magnitude'][bright_stars]
marker_size = (0.5 + limiting_magnitude - magnitude) ** 2.0

xy1 = stars[['ra_degrees', 'dec_degrees']].loc[edges_star1].values
xy2 = stars[['ra_degrees', 'dec_degrees']].loc[edges_star2].values
lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)


ras = np.array(ras)
decs = np.array(decs)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)
ax.set_xlabel('Longitude')
ax.set_ylabel('Altitude')
ax.set_title('Observation Chart')

unique_longitudes = np.unique(ras)
colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_longitudes)))

for lon, color in zip(unique_longitudes, colors):
    mask = ras == lon
    ax.scatter(ras[mask], decs[mask], label=f'Longitude {lon}', color=color)

ax.add_collection(LineCollection(lines_xy, colors='#00f2'))
ax.scatter(stars['ra_degrees'][bright_stars], stars['dec_degrees'][bright_stars],
           s=marker_size, color='k')

ax.grid(True)

for i, name in enumerate(celestial_bodies):
    ax.annotate(name.split()[0], (ras[i], decs[i]), textcoords='offset points', xytext=(5, 5), ha='center')

#fig.savefig('astro.png', bbox_inches='tight')
plt.show()