import numpy as np
from environment_mod.geo_lat2geod_lat import geo_lat2geod_lat
def geo2geod(lat_geo_phi, lon_geo_lmd, alt_geo):
    lat_geod_phi = geo_lat2geod_lat(lat_geo_phi)  # degrees
    lon_geod_lmd = lon_geo_lmd  # degrees
    alt_geod = alt_geo  # km
    return lat_geod_phi, lon_geod_lmd, alt_geod