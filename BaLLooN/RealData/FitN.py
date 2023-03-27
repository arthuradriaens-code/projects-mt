import NuRadioReco.detector.detector
import numpy as np
import radiotools
from NuRadioReco.utilities import units
import NuRadioReco.modules.io.rno_g.readRNOGData
from NuRadioReco.modules import channelAddCableDelay
from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
from NuRadioMC.utilities import medium
from NuRadioMC.SignalProp import radioproparaytracing
import matplotlib.pyplot as plt
import os
from coordinate_system import CoordinateSystem
import datetime
import folium
from dateutil.tz import tzutc
import copy
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

station_id = 23
channel_a_id = 5
channel_b_id = 6

#-------------------------------------------------------------------------------#
#                               functions                                       #
#-------------------------------------------------------------------------------#

def Sine(t,A,T,offset):
    return A*np.sin(0.403*2*np.pi*t + T) + offset
def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degreetorad(deg):
    return deg*np.pi/180
def degrees(SomethingInRads):
    return SomethingInRads*180/np.pi
def delta_taccent(theta,deltaz,n):
    v = c/n
    return ((np.cos(theta)*deltaz)/v)*(10**9)

#-------------------------------------------------------------------------------#
#                               constants                                       #
#-------------------------------------------------------------------------------#
c = 299792458 #(m/s)
ice = medium.greenland_simple()
indexofrefractionrange = np.linspace(1.27,2.5,10000)
n_icesurface = ice.get_index_of_refraction(np.array([0,0,-0.00001]))

#-------------------------------------------------------------------------------#
#                               spatial data                                    #
#-------------------------------------------------------------------------------#
gpx_file = "/home/arthur/Documents/thesis/data/sonde/gpx/SMT_20220829_111459.gpx"
StationNumber = 23
GivenTime = "2022/08/29/11/18/32"
i = 0
track = 0
for track in read_gpx_file(gpx_file):
    continue
segment = track['segments'][0]

masked_segment = {key:[] for key in segment.keys()}
masked_segment['elevation-up'] = segment['elevation-up']
masked_segment['elevation-down'] = segment['elevation-down']

for i, time in enumerate(segment['time']):
    for key in list(segment.keys())[:-2]:
        masked_segment[key].append(segment[key][i])

for key in list(segment.keys())[:-2]:
    if key in ['time','latlon']: continue
    masked_segment[key] = np.array(masked_segment[key])

GivenTime = GivenTime.split('/')
GivenTime = datetime.datetime(int(GivenTime[0]), int(GivenTime[1]), int(GivenTime[2]),int(GivenTime[3]),int(GivenTime[4]),int(GivenTime[5]),tzinfo=tzutc())
print("looking at the event recorded at:")
print(GivenTime)
coor = CoordinateSystem()
for t,Time in enumerate(segment['time']):
    if (Time - GivenTime).total_seconds() < 1: 
        BalloonPosition = coor.geodetic_to_enu(segment['lat'][t],segment['lon'][t],segment['elevation'][t])
print("balloon position:")
print(BalloonPosition)

# data reader
data_reader = NuRadioReco.modules.io.rno_g.readRNOGData.readRNOGData()
data_reader.begin("/home/arthur/Documents/thesis/data/interesting/station23/run691/combined.root")


# get event
for event in data_reader.run():
    if event.get_id()==489:
        print("found you")
        break

# get detector information at the specified event time
det = NuRadioReco.detector.detector.Detector(json_filename="/home/arthur/Documents/thesis/programs/analysis-tools/rnog_analysis_tools/detector_json/RNO_season_2022.json", 
                                             antenna_by_depth=False)
station = event.get_station(station_id)
det.update(station.get_station_time())
channel_a = station.get_channel(channel_a_id)
channel_b = station.get_channel(channel_b_id)
cable_delay_a = det.get_cable_delay(station_id,channel_a_id)
cable_delay_b = det.get_cable_delay(station_id,channel_b_id)

#detector
AddCableDelay = channelAddCableDelay.channelAddCableDelay()
# channels:
locallocationstation = np.array(det.get_absolute_position(station_id))
print(locallocationstation)

# Ok now we have our balloon position and the detector position, 
# We'll set the BalloonPosition relative to the detector:

Balloon = np.array([BalloonPosition[0]-locallocationstation[0],BalloonPosition[1] - locallocationstation[1],BalloonPosition[2]] - locallocationstation[2])

#And let's rotate to only have an x component, no y:
r = np.sqrt(Balloon[0]**2 + Balloon[1]**2)
Balloon[1] = 0
Balloon[0] = r


#-------------------------------------------------------------------------------#
#                           Get observed time difference                        #
#-------------------------------------------------------------------------------#
# NOTE: This doesn't work yet

Detectors = []
print(det.get_relative_position(station_id,channel_a_id))
Detectors.append(np.array(det.get_relative_position(station_id,channel_a_id)) * units.m)
Detectors.append(np.array(det.get_relative_position(station_id,channel_b_id)) * units.m)
print(det.get_relative_position(station_id,channel_b_id))

# channel a fit
channel_a_voltages = channel_a.get_trace()
channel_a_spectrum = channel_a.get_frequency_spectrum()
channel_a_frequencies = channel_a.get_frequencies()
plt.axvline(x=0.403,linestyle='dashed')
plt.plot(channel_a_frequencies,channel_a_spectrum)
plt.show()
channel_a_times = channel_a.get_times()
plt.plot(channel_a_times,channel_a_voltages)
fitted_a, cov_fitted_a = curve_fit(Sine, channel_a_times,channel_a_voltages)
print(fitted_a)
plt.plot(channel_a_times,Sine(channel_a_times,*fitted_a))
plt.show()

# channel b fit
channel_b_spectrum = channel_b.get_frequency_spectrum()
channel_b_frequencies = channel_b.get_frequencies()
plt.axvline(x=0.403,linestyle='dashed')
plt.plot(channel_b_frequencies,channel_b_spectrum)
plt.show()
channel_b_voltages = channel_b.get_trace()
channel_b_times = channel_b.get_times()
plt.plot(channel_b_times,channel_b_voltages)
fitted_b, cov_fitted_b = curve_fit(Sine, channel_b_times,channel_b_voltages)
plt.plot(channel_b_times,Sine(channel_b_times,*fitted_b))
plt.show()

difference = np.abs(fitted_b[1] - fitted_a[1] - cable_delay_a  + cable_delay_b)
# time is, for some reason in hunderds of nanoseconds
print(difference)

NumberOfDetectors = len(Detectors)
delta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))
delta_t[0][1] = difference
deltaz = np.linalg.norm(Detectors[0]-Detectors[1])
if difference < deltaz*0.3:
    print("The signal appears to be moving faster than light...")
if difference > deltaz*0.6:
    print("The signal moves WAY to slow")
#-------------------------------------------------------------------------------#
#                                   Fit n                                       #
#-------------------------------------------------------------------------------#
differences = np.zeros(len(indexofrefractionrange))

b_ballon = -67.704
a_ballon = (Balloon[2]-b_ballon)/Balloon[0]

for number,n in enumerate(indexofrefractionrange):
    thetas = np.linspace(0,np.pi/2,1000)
    delta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    correlation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    normedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    summedcorrelation = np.zeros(1000)

    for i in range(NumberOfDetectors):
        for j in range(NumberOfDetectors):
            if i < j:
                deltaz = np.linalg.norm(Detectors[i]-Detectors[j])
                delta_taccenten[i][j] = delta_taccent(thetas,np.abs(deltaz),n)
                correlation[i][j] = np.abs(delta_t[i][j] - delta_taccenten[i][j])
                normedcorrelation[i][j] = correlation[i][j]/np.trapz(correlation[i][j],thetas)

                summedcorrelation += normedcorrelation[i][j]

    angle_index = np.where(summedcorrelation == summedcorrelation.min())
    angle = thetas[angle_index] #zenith ofc

    a_planewave = np.tan(np.pi/2-angle)
    b_planewave = -67.704

    angle_snell = np.arcsin(np.sin(angle)*n_icesurface)
    a_snell = np.tan(np.pi/2-angle_snell)
    b_snell = -1*a_snell*(-1*b_planewave/a_planewave)
    XopBallonHoogte = (Balloon[2] - b_snell)/a_snell
    verschil = XopBallonHoogte - Balloon[0]

    differences[number] = np.abs(verschil)

print("minimal difference between timing:")
print(differences.min())
n_index = np.where(differences == differences.min())
n_fit = indexofrefractionrange[n_index]
if len(n_fit) > 1:
    print("undetermined")
    print(n_fit)
    n_fit = n_fit[0]
print("index of refraction from fit: {}".format(n_fit))


