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
    ra_deg = ra.dstr(warn=False).split('deg')[0].strip()
    dec_deg = dec.dstr(warn=False).split('deg')[0].strip()
    print(f"\n{body}:")
    print(ra_deg)
    print(dec_deg)
    ras.append(int(ra_deg))
    decs.append(int(dec_deg))

"""
Stars and Constelations
"""

# Choose the constellations to plot
names = ["Cru", "Ari", "Tau", "Gem", "Cnc", "Leo", "Vir", "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

with load.open(hipparcos.URL) as f:
    stars = hipparcos.load_dataframe(f)

url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
       '/skycultures/modern_st/constellationship.fab')

with load.open(url) as f:
    constellations = stellarium.parse_constellations(f)

edges = [edge for name, edges in constellations for edge in edges if name in names]
edges_star1 = [star1 for star1, star2 in edges]
edges_star2 = [star2 for star1, star2 in edges]

alt, az, _ = observer_on_earth.at(t).altaz()
center = observer.at(t).from_altaz(alt_degrees=alt.degrees, az_degrees=az.degrees)
projection = build_stereographic_projection(center)

star_positions = observer_on_earth.at(t).observe(Star.from_dataframe(stars))

limiting_magnitude = 4.40

bright_stars = (stars.magnitude <= limiting_magnitude)

magnitude = stars['magnitude'][bright_stars]
marker_size = (0.5 + limiting_magnitude - magnitude) ** 2.0

# Convert RA and Dec to radians
def convert_ra(values):
    values = list(values)
    return [(value + 180) % 360 - 180 for value in values]

ra_rad = np.deg2rad(convert_ra(ras))
dec_rad = np.deg2rad(decs)

stars['ra_degrees'] = convert_ra(stars['ra_degrees'])

stars_ra_rad = np.deg2rad(stars['ra_degrees'][bright_stars])
stars_dec_rad = np.deg2rad(stars['dec_degrees'][bright_stars])

xy1 = np.deg2rad(stars[['ra_degrees', 'dec_degrees']].loc[edges_star1]).values
xy2 = np.deg2rad(stars[['ra_degrees', 'dec_degrees']].loc[edges_star2]).values
lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111, projection='aitoff')

ax.scatter(stars_ra_rad, stars_dec_rad,
           s=marker_size, color='k')

unique_longitudes = np.unique(ra_rad)
colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_longitudes)))

for lon, color in zip(unique_longitudes, colors):
    mask = ra_rad == lon
    ax.scatter(ra_rad[mask], dec_rad[mask], label=f'Longitude {lon}', color=color)
    
for i, name in enumerate(celestial_bodies):
    ax.annotate(name.split()[0], (ra_rad[i], dec_rad[i]), textcoords='offset points', xytext=(5, 5), ha='center')


ax.add_collection(LineCollection(lines_xy, colors='#00f2'))

ax.grid(True)
ax.set_title('Galactic Lat-Long Plot')
ax.set_xlabel('Galactic Longitude (radians)')
ax.set_ylabel('Galactic Latitude (radians)')

fig.savefig('sky_observation.png', bbox_inches='tight')
plt.show()