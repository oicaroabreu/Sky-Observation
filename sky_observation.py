import matplotlib.pyplot as plt
from skyfield.api import load
from skyfield.api import N,S,E,W, wgs84
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
observer = earth + wgs84.latlon(latitude * N, longitude * W)

altitudes = []
azimuths = []

for body in celestial_bodies:
    planet = planets[body]
    planet_pos = observer.at(t).observe(planet)
    app = planet_pos.apparent()
    alt, az, _ = app.altaz()
    print(f"\n{body}:")
    print(alt.dstr())
    print(az.dstr())
    altitudes.append(alt.degrees)
    azimuths.append(az.degrees)

altitudes = np.array(altitudes)
azimuths = np.array(azimuths)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(np.min(azimuths) - 20, np.max(azimuths)+ 20)
ax.set_ylim(np.min(altitudes), np.max(altitudes))
ax.set_xlabel('Longitude')
ax.set_ylabel('Altitude')
ax.set_title('Observation Chart')

unique_longitudes = np.unique(azimuths)
colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_longitudes)))

for lon, color in zip(unique_longitudes, colors):
    mask = azimuths == lon
    ax.scatter(azimuths[mask], altitudes[mask], label=f'Longitude {lon}', color=color)

ax.grid(True)

for i, name in enumerate(celestial_bodies):
    ax.annotate(name, (azimuths[i], altitudes[i]), textcoords='offset points', xytext=(5, 5), ha='center')

plt.show()