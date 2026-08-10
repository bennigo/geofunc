"""Pin the axis order of :func:`geofunc.geo.xyzell`.

The docstring claimed ``(lon, lat, height)`` for years while the function
actually returned ``(lat, lon, height)``. That was true under pyproj 1.x, where
``lonlat`` was ``proj.Proj(init="EPSG:4326")`` and traditional x/y order
applied; the pyproj 1→2 migration to ``CRS("EPSG:4326")`` silently flipped it
to EPSG:4326's authority order (lat, lon) and nothing caught it.

These tests exist so a future pyproj upgrade — or someone adding
``always_xy=True`` — cannot flip it back unnoticed. Iceland is a good fixture
precisely because latitude (~64 N) and longitude (~-22 E) cannot be confused
for one another.
"""

import numpy as np
import pytest

from geofunc.geo import xyzell

# REYK, from the IMO a priori coordinate file (ITRF2008 / GRS80).
REYK_XYZ = [2587383.91122, -1043033.58722, 5716564.1974]
REYK_LAT = 64.138788
REYK_LON = -21.955490

# ELDC (Eldvörp), the station that exposed the bug in production.
ELDC_XYZ = [2602396.65469, -1079316.76219, 5703115.26511]
ELDC_LAT = 63.863143
ELDC_LON = -22.525718


@pytest.mark.parametrize(
    "xyz,lat,lon",
    [(REYK_XYZ, REYK_LAT, REYK_LON), (ELDC_XYZ, ELDC_LAT, ELDC_LON)],
)
def test_returns_lat_lon_height_in_that_order(xyz, lat, lon):
    out = xyzell(xyz, radians=False)

    assert out[0] == pytest.approx(lat, abs=1e-5), "element 0 must be LATITUDE"
    assert out[1] == pytest.approx(lon, abs=1e-5), "element 1 must be LONGITUDE"


def test_latitude_and_longitude_are_not_swapped():
    """A blunt guard: Iceland's lat and lon are not interchangeable."""
    lat, lon, _ = xyzell(REYK_XYZ, radians=False)

    assert 60.0 < lat < 70.0, f"element 0 ({lat}) is not an Icelandic latitude"
    assert -30.0 < lon < -10.0, f"element 1 ({lon}) is not an Icelandic longitude"


def test_height_is_ellipsoidal_not_orthometric():
    """Iceland's geoid separation is ~60-70 m; an orthometric height would
    differ from the ellipsoidal one by roughly that much."""
    _, _, h = xyzell(REYK_XYZ, radians=False)

    assert h == pytest.approx(93.028, abs=0.01)


def test_radians_is_the_default():
    rad = xyzell(REYK_XYZ)
    deg = xyzell(REYK_XYZ, radians=False)

    assert np.degrees(rad[0]) == pytest.approx(deg[0], abs=1e-9)
    assert np.degrees(rad[1]) == pytest.approx(deg[1], abs=1e-9)
    # Height is metres in both cases, untouched by the radians flag.
    assert rad[2] == pytest.approx(deg[2], abs=1e-9)
