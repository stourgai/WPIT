def geo_lat2geod_lat(phi):
    a = 6378137  # meter semi major axis of earth
    f = 1 / 298.257  # flattening
    b = a - f * a  # semi minor axis
    e = ((a ** 2 - b ** 2) ** (1 / 2)) / a
    phi_rad = np.deg2rad(phi)
    geod_lat = np.arctan(np.tan(phi_rad) / (1 - e ** 2))
    geod_lat = np.rad2deg(geod_lat)
    return geod_lat  # in degrees
