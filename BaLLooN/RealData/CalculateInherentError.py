# Given a balloon file, a detector and the time, this will calculate the
# systematic error that's to be expected


from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
import os
from coordinate_system import CoordinateSystem
import datetime
import folium
from dateutil.tz import tzutc
import copy
import numpy as np
import csv


#gpx_file_name = 'SMT_20220630_112143.gpx'  #Jess
#gpx_file_name = 'SMT_20220701_231934.gpx' #Bob
#gpx_file_name = 'SMT_20220715_231621.gpx' #Bob
print("give the filename, e.g /data/sonde/gpx/SMT_20220701_231934.gpx")
gpx_file = input()
i = 0
track = 0
for track in read_gpx_file(gpx_file):
    continue

segment = track['segments'][0]

print("What station? E.g 21")
StationNumber = input()
stations = {21:[72.5874063909459,-38.4660301212611], 12:[72.6000868058195,-38.4962265332872],11:[72.5892267215905,-38.5022988244688], 13:[72.6109470001738,-38.4901465440588], 22:[72.598265271346,-38.4599355034766], 23:[72.6091242603966,-38.4538331609837] ,24:[72.6199833575357,-38.4477230792255]}
StationCoordinate = stations[int(StationNumber)]
coor = CoordinateSystem()
locallocationstation = coor.geodetic_to_enu(StationCoordinate[0],StationCoordinate[1],3251.9489147560234)

masked_segment = {key:[] for key in segment.keys()}
masked_segment['elevation-up'] = segment['elevation-up']
masked_segment['elevation-down'] = segment['elevation-down']

for i, time in enumerate(segment['time']):
    for key in list(segment.keys())[:-2]:
        masked_segment[key].append(segment[key][i])

for key in list(segment.keys())[:-2]:
    if key in ['time','latlon']: continue
    masked_segment[key] = np.array(masked_segment[key])

print("at what time? E.g 2022/07/24/23/19/27")
GivenTime = input()
GivenTime = GivenTime.split('/')
GivenTime = datetime.datetime(int(GivenTime[0]), int(GivenTime[1]), int(GivenTime[2]),int(GivenTime[3]),int(GivenTime[4]),int(GivenTime[5]),tzinfo=tzutc())
print(GivenTime)
print(segment['time'][5])
i = np.where((segment['time'] - GivenTime).total_seconds() < 1)
BalloonPosition = coor.geodetic_to_enu(segment['lat'][i],segment['lon'][i],segment['elevation'][i])
print("balloon height:")
print(BalloonPosition[2])

# Ok now we have our balloon position and the detector position, let's use the simulation 
# part to figure out what the systematic error is.

Detectorx = locallocationstation[0]
Detectory = locallocationstation[1]
Detectors = np.zeros((4,3))
Detectors[0] = np.array([Detectorx, Detectory, -97.]) * units.m
Detectors[1] = np.array([Detectorx, Detectory, -96.]) * units.m
Detectors[2] = np.array([Detectorx, Detectory, -95.]) * units.m
Detectors[3] = np.array([Detectorx, Detectory, -94.]) * units.m
# Middle of channels 0-3, making sure about the depth
traveltimes = []
paths = []
times = []
distances = []
Balloon = BalloonPosition
prop = radioproparaytracing.radiopropa_ray_tracing(ice, attenuation_model='GL1',config=configh)
for detector in Detectors:
    start_point = Balloon
    final_point = detector
    prop.set_start_and_end_point(start_point, final_point)
    prop.find_solutions()
    SolNumber = prop.get_number_of_solutions()
    for Sol in range(SolNumber):
        paths.append(prop.get_path(Sol))
        times.append(prop.get_travel_time(Sol))
        x = np.linspace(paths[-1][0,0],start_point[0],1000)
        xlen = x[-1]-x[0]
        z = np.linspace(paths[-1][0,2],start_point[2],1000)
        zlen = z[-1]-z[0]
        diagonallen = np.sqrt(xlen*xlen + zlen*zlen) #(m)
        traveltime = times[-1]/units.ns + (diagonallen/c)*(10**9) #ns
        traveltimes.append(traveltime)

def delta_taccent(theta,deltaz,n):
    v = c/n
    return ((np.cos(theta)*deltaz)/v)*(10**9)

differences = np.zeros(len(indexofrefractionrange))

for number,n in enumerate(indexofrefractionrange):
    thetas = np.linspace(0,np.pi/2,1000)
    NumberOfDetectors = len(Detectors)
    delta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))
    delta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    correlation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    normedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    summedcorrelation = np.zeros(1000)

    for i in range(NumberOfDetectors):
        for j in range(NumberOfDetectors):
            if i < j:
                delta_t[i][j] = traveltimes[i] - traveltimes[j]
                deltaz = np.linalg.norm(Detectors[i]-Detectors[j])
                delta_taccenten[i][j] = delta_taccent(thetas,np.abs(deltaz),n)
                correlation[i][j] = np.abs(delta_t[i][j] - delta_taccenten[i][j])
                normedcorrelation[i][j] = correlation[i][j]/np.trapz(correlation[i][j],thetas)

                summedcorrelation += normedcorrelation[i][j]

    angle_index = np.where(summedcorrelation == summedcorrelation.min())
    angle = thetas[angle_index] #zenith ofc

    x_relative = Balloon[0] - Detectors[0][0]
    y_relative = Balloon[1] - Detectors[0][1]
    z_relative = Balloon[2] + 95.5
    direct_angle = np.arccos(z_relative/(x_relative**2 + y_relative**2 + z_relative**2))
    differences[number] = np.abs(direct_angle - angle)

n_index = np.where(differences == differences.min())
n_fit = indexofrefractionrange[n_index]
if len(n_fit) > 1:
    print("undetermined")
    print(n_fit)
    n_fit = n_fit[0]
print("index of refraction from fit: {}".format(n_fit))

position = np.array([0,0,-95.5])
direct_angle = np.pi/2 - np.arctan(a_ballon)

n_actual = ice.get_index_of_refraction(position)
print("actual index of refraction from model: {}".format(n_actual))

print("relative accuracy: {}%".format(100*(n_fit - n_actual)/n_actual))
