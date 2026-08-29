"""geofunc — pure-Python geodesy primitives for the GPS library ecosystem.

Tier 0 foundation library. Coordinate reference-frame conversions live in
:mod:`geofunc.coords`; geodesy math (ellipsoidal/Cartesian conversions,
plate motion) in :mod:`geofunc.geo`.
"""

__path__ = __import__("pkgutil").extend_path(__path__, __name__)

from .coords import CoordFrame, FRAMES, frame_names, get_frame, transform

__all__ = [
    "CoordFrame",
    "FRAMES",
    "frame_names",
    "get_frame",
    "transform",
]


def main() -> None:
    print("Hello from geofunc!")
