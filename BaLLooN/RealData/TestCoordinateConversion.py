from coordinate_system import CoordinateSystem

coor = CoordinateSystem()
print(coor.geodetic_to_enu(72.6000868058195,-38.4962265332872,3251.9489147560234))
print(coor.geodetic_to_enu(72.6000868058195,-38.4962265332872,3251.9489147560234)[0])
