# geofunc

Pure-Python geodesy primitives for the GPS library ecosystem.
**Tier 0** — foundation library, zero internal dependencies.

## Status

Stable. ~649 LOC across two modules. No active feature work; consumed by `geo_dataread`.

## Layout

```
geofunc/
├── src/geofunc/
│   ├── __init__.py       # public API: re-exports coords (transform, FRAMES, …)
│   ├── geo.py            # core geodesy: ellipsoidal/Cartesian conversions, plate motion
│   ├── geofunc.py        # higher-level helpers: station coords, plate-velocity wrappers
│   └── coords.py         # reference-frame registry + transform() facade (ITRF/WGS84/ISN)
├── tests/                # test_geo.py, test_coords.py + ad-hoc scripts
├── perlscripts/          # legacy Perl utilities (reference only)
└── pyproject.toml
```

## Key Functions

| Module | Function | Purpose |
|--------|----------|---------|
| `geo.py` | `xyzell` | ECEF (X,Y,Z) → ellipsoidal (lat, lon, h) |
| `geo.py` | `platePoles` | Tectonic plate Euler poles |
| `geo.py` | `rotpole` | Rotate pole between reference frames |
| `geo.py` | `greatcircd` | Great-circle distance |
| `geofunc.py` | `convllh` | Coordinate conversion driver |
| `geofunc.py` | `plateVelo` | Plate velocity at a point |
| `geofunc.py` | `getStationCoordinates` | Look up station coordinates |
| `coords.py` | `transform` | Convert (x,y,z) between named frames (wgs84, itrf08/14/20, isn93/2004/2016) |
| `coords.py` | `FRAMES` / `get_frame` | Data-driven frame registry (add a frame = one data row) |

## Dependencies

- **In** (this package needs): `pyproj`, `numpy`, `scipy` (external only).
  `geofunc.py` additionally imports `gps_parser` (config-driven helpers) — the
  `coords.py` conversion core is deliberately pure (pyproj + numpy only).
- **Out** (used by): `geo_dataread` (coordinate transforms, plate-velocity calcs).

### Coordinate frames (`coords.py`)

Static frame-to-frame only (no epoch/plate-motion — that lands with the data
analyses). Registered frames: `wgs84` (lon/lat/h), `itrf2008`/`itrf2014`/
`itrf2020` (ECEF X/Y/Z), `isn93`/`isn2004`/`isn2016` (Lambert E/N). Axis order
is always x-first (`always_xy=True`), unlike `geo.xyzell`'s (lat, lon, h).
Adding a frame is a data row in `coords._FRAME_SPECS`. `tos search` will consume
this as computed coordinate columns (separate work).

## Console Scripts

```bash
geofunc       # entry: geofunc:main
```

## Cross-References

- `../CLAUDE.md` — ecosystem overview + dependency graph
- Vault hub: `/home/bgo/notes/bgovault/2.Areas/VI_GPS_Library/1776347706-gps-library-ecosystem-hub.md`

---

*Last reviewed: 2026-08-29*
