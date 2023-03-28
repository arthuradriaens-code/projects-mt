# Given a balloon file, a detector and the time, this will calculate the
# systematic error that's to be expected


from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
from NuRadioMC.utilities import medium
from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import os
from coordinate_system import CoordinateSystem
import datetime
import folium
from dateutil.tz import tzutc
import copy
import numpy as np
import csv

def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degreetorad(deg):
    return deg*np.pi/180
def degrees(SomethingInRads):
    return SomethingInRads*180/np.pi

c = 299792458 #(m/s)
ice = medium.greenland_simple()
indexofrefractionrange = np.linspace(1.4,2,10000)
n_icesurface = ice.get_index_of_refraction(np.array([0,0,-0.00001]))

configh = dict()
configh['propagation'] = dict(
    attenuate_ice = True,
    focusing_limit = 2,
    focusing = False,
    radiopropa = dict(
        mode = 'hybrid minimizing',
        iter_steps_channel = [45., 2., .5, .05,0.005], #unit is meter
        iter_steps_zenith = [.7, .05, .005, .001,0.0001], #unit is degree
        auto_step_size = False,
        max_traj_length = 10000) #unit is meter
)
configh['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)


#gpx_file_name = 'SMT_20220630_112143.gpx' #Jess
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
locallocationstation = coor.geodetic_to_enu(StationCoordinate[0],StationCoordinate[1])
print(locallocationstation)

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
for t,Time in enumerate(segment['time']):
    if (Time - GivenTime).total_seconds() < 1: 
        BalloonPosition = coor.geodetic_to_enu(segment['lat'][t],segment['lon'][t],segment['elevation'][t])
print("balloon position:")
print(BalloonPosition)

# Ok now we have our balloon position and the detector position, let's use the simulation 
# part to figure out what the systematic error is.
# We'll set the BalloonPosition relative to the detector:

Balloon = np.array([BalloonPosition[0]-locallocationstation[0],BalloonPosition[1] - locallocationstation[1],BalloonPosition[2]] - locallocationstation[2])

#And let's rotate to only have an x component, no y:
r = np.sqrt(Balloon[0]**2 + Balloon[1]**2)
Balloon[1] = 0
Balloon[0] = r

# channels:
Detectorx = 0
Detectory = 0
#Detectors = np.zeros((4,3))
Detectors = np.zeros((2,3))
Detectors[0] = np.array([Detectorx, Detectory, -37.719]) * units.m
Detectors[1] = np.array([Detectorx, Detectory, -57.709]) * units.m
#Detectors[2] = np.array([Detectorx, Detectory, -95.]) * units.m
#Detectors[3] = np.array([Detectorx, Detectory, -94.]) * units.m

traveltimes = []
paths = []
times = []
distances = []


prop = radioproparaytracing.radiopropa_ray_tracing(ice, attenuation_model='GL1',config=configh)
for detector in Detectors:
    start_point = Balloon
    print(start_point)
    final_point = detector
    print(final_point)
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

        plt.plot(x,z,color='orange')
        plt.plot(paths[-1][:,0],paths[-1][:,2],label="travel time = {0:.2f} nanoseconds".format(traveltime) ,color="orange")
def delta_taccent(theta,deltaz,n):
    v = c/n
    return ((np.cos(theta)*deltaz)/v)*(10**9)

differences = np.zeros(len(indexofrefractionrange))

b_ballon = -47.7
a_ballon = (Balloon[2]-b_ballon)/Balloon[0]

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

    a_planewave = np.tan(np.pi/2-angle)
    b_planewave = -47.7

    angle_snell = np.arcsin(np.sin(angle)*n_icesurface)
    a_snell = np.tan(np.pi/2-angle_snell)
    b_snell = -1*a_snell*(-1*b_planewave/a_planewave)
    XopBallonHoogte = (Balloon[2] - b_snell)/a_snell
    verschil = XopBallonHoogte - Balloon[0]

    differences[number] = np.abs(verschil)

print(differences.min())
n_index = np.where(differences == differences.min())
n_fit = indexofrefractionrange[n_index]
if len(n_fit) > 1:
    print("undetermined")
    print(n_fit)
    n_fit = n_fit[0]
print("index of refraction from fit: {}".format(n_fit))

position = np.array([0,0,-47.7])

n_actual = ice.get_index_of_refraction(position)
print("actual index of refraction from model: {}".format(n_actual))

#plotting
n = n_fit
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
a_planewave = np.tan(np.pi/2-angle)
b_planewave = -47.7

angle_snell = np.arcsin(np.sin(angle)*n_icesurface)
a_snell = np.tan(np.pi/2-angle_snell)
b_snell = -1*a_snell*(-1*b_planewave/a_planewave)

b_ballon = -47.7
a_ballon = (Balloon[2]-b_ballon)/Balloon[0]
x = np.linspace(0,Balloon[0],1000)
plt.plot(x,a_ballon*x + b_ballon,color="blue")
a_refracted = np.tan(np.pi/2-angle)
b_refracted = -47.7

y = a_refracted*x + b_refracted
indexwhensurface = np.where(y > 0)[0][0]
x = x[0:indexwhensurface]
y = y[0:indexwhensurface]
plt.plot(x,y,color="red",label="Plane Wave Reconstructed Path")
x_snell = np.linspace(x[-1],Balloon[0])
y = a_snell*x_snell + b_snell
plt.plot(x_snell,y,color="red")

plt.ylabel("vertical distance (m)")
plt.xlabel("horizontal distance (m)")
plt.title("Greenland simple trajectory with GL1 attenuation\n solved with hybrid ray tracer")
plt.ylim(-100,Balloon[0]+5)
plt.xlim(-100,Balloon[2]+5)
plt.legend()
plt.show()

print("relative accuracy: {}%".format(100*(n_fit - n_actual)/n_actual))
