import sys
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
from scipy.signal import find_peaks

station_id = 21
event_id = 117
channel_ids = [1,2]
traveltimes = np.zeros(len(channel_ids))
gpx_file = "/mnt/usb/sonde/gpx/SMT_20220726_111605.gpx"
GivenTime = "2022/07/26/11/18/40"
rootfile = "/mnt/usb/RNO-G-DATA/station21/run1441/combined.root"
#phi2s = [-786.91,-894.548,-895,-782]
phi2s = [-636.55,-630.964]
#phi2s = [-3920.707,-3918.697]
#phi2s = [-727.8,-658]
#phi2s = [-786.91,-782]

#-------------------------------------------------------------------------------#
#                               Colors                                          #
#-------------------------------------------------------------------------------#

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#-------------------------------------------------------------------------------#
#                               functions                                       #
#-------------------------------------------------------------------------------#

def AM(t,A,m,frel,phi1,fm,phi2):
    return A*np.sin(0.403125*frel*units.ns*2*np.pi*t + phi1)*(m*np.cos(fm*2*np.pi*t + phi2) + 1)
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
def hl_envelopes_idx(s, dmin=1, dmax=1, split=False):
    """
    Input :
    s: 1d-array, data signal from which to extract high and low envelopes
    dmin, dmax: int, optional, size of chunks, use this if the size of the input signal is too big
    split: bool, optional, if True, split the signal in half along its mean, might help to generate the envelope in some cases
    Output :
    lmin,lmax : high/low envelope idx of input signal s
    """

    # locals min      
    lmin = (np.diff(np.sign(np.diff(s))) > 0).nonzero()[0] + 1 
    # locals max
    lmax = (np.diff(np.sign(np.diff(s))) < 0).nonzero()[0] + 1 
    
    if split:
        # s_mid is zero if s centered around x-axis or more generally mean of signal
        s_mid = np.mean(s) 
        # pre-sorting of locals min based on relative position with respect to s_mid 
        lmin = lmin[s[lmin]<s_mid]
        # pre-sorting of local max based on relative position with respect to s_mid 
        lmax = lmax[s[lmax]>s_mid]

    # global min of dmin-chunks of locals min 
    lmin = lmin[[i+np.argmin(s[lmin[i:i+dmin]]) for i in range(0,len(lmin),dmin)]]
    # global max of dmax-chunks of locals max 
    lmax = lmax[[i+np.argmax(s[lmax[i:i+dmax]]) for i in range(0,len(lmax),dmax)]]
    
    return lmin,lmax
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
#-------------------------------------------------------------------------------#
#                               constants                                       #
#-------------------------------------------------------------------------------#
c = 299792458 #(m/s)
ice = medium.greenland_simple()
indexofrefractionrange = np.linspace(1.3,1.8,50000)
n_icesurface = ice.get_index_of_refraction(np.array([0,0,-0.00001]))

#-------------------------------------------------------------------------------#
#                               spatial data                                    #
#-------------------------------------------------------------------------------#
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

# data reader
data_reader = NuRadioReco.modules.io.rno_g.readRNOGData.readRNOGData()
data_reader.begin(rootfile)


# get event
for event in data_reader.run():
    if event.get_id()==event_id:
        print("found you")
        break

# get detector information at the specified event time
det = NuRadioReco.detector.detector.Detector(json_filename="/home/arthur/Universiteit/master-proef/analysis-tools/rnog_analysis_tools/detector_json/RNO_season_2022.json", 
                                             antenna_by_depth=False)
station = event.get_station(station_id)
det.update(station.get_station_time())
passband = [0.40312499,0.40312501]


locallocationstation = np.array(det.get_absolute_position(station_id))

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
# NOTE: This doesn't seem to work yet


Detectors = []
#approxtimes = np.array([9774,9779,9785,9790])
#omegas = np.array([0.16,0.152,0.151,0.162])
#delays = np.array([700.65,705,710.3,714.96])
#phis = (approxtimes - delays)*omegas
#phis = np.zeros(len(channel_ids))
#phis[0] = -611.02040816
#phis[1] = -615.30612245
#phis[3] = 1471
#print(phis)

for i,channel_id in enumerate(channel_ids):
    channel = station.get_channel(channel_id)
    cable_delay = det.get_cable_delay(station_id,channel_id)
    Detectors.append(np.array(det.get_relative_position(station_id,channel_id)) * units.m)
    channel_voltages = channel.get_filtered_trace(passband, 'butter', 6) 
    channel_spectrum = channel.get_frequency_spectrum()
    channel_frequencies = channel.get_frequencies()
    plt.axvline(x=0.403125,linestyle='dashed',color="grey")
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Counts")
    plt.title("Frequency spectrum of channel {}".format(channel_id))
    plt.plot(channel_frequencies,channel_spectrum)
    plt.show()
    channel_times = channel.get_times()
    plt.plot(channel_times,channel_voltages,label="measured data")
    lmin,lmax = hl_envelopes_idx(channel_voltages)
    plt.plot(channel_times[lmin],channel_voltages[lmin],label="low envelope")
    plt.plot(channel_times[lmax],channel_voltages[lmax],label="high envelope")
    E_1=np.max(channel_voltages[lmax])
    E_2=np.abs(np.max(channel_voltages[lmin]))
    m = (E_1/E_2 - 1)/(E_1/E_2 + 1)
    A = E_1/(1+m)
    listofmaxes,_ = find_peaks(channel_voltages[lmax])
    fm = 1/((channel_times[lmax][listofmaxes[1]] - channel_times[lmax][listofmaxes[0]]))
    frel=1
    phi1 = 0
    phi2 = phi2s[i]
    fitted, cov_fitted = curve_fit(AM,channel_times,channel_voltages,p0=(A,m,frel,phi1,fm,phi2))
    phase_2 = fitted[-1]
    T2 = np.abs(phase_2/(2*np.pi*fitted[-2]))
    traveltimes[i] = np.abs(T2 - cable_delay)
    print("\n")
    print("-----------------------------")
    print("Channel {}".format(channel_id))
    print("-----------------------------")
    print("phi2 used: {}".format(phase_2))
    print("omega used: {}".format(2*np.pi*fitted[-2]))
    print("cable delay: {}".format(cable_delay))
    print("T2: {}".format(T2))
    print("traveltime: {}".format(traveltimes[i]))
    print("-----------------------------")
    print("\n")
    phase_1 = fitted[-3]
    T1 = phase_1/(2*np.pi*0.403125)
    plt.xlabel("time (nanoseconds)")
    plt.ylabel("Voltage")
    plt.title("Voltage i.f.o time for channel {}".format(channel_id))
    plt.plot(channel_times,AM(channel_times,*fitted),label="fitted sine",linestyle="dotted")
    plt.legend()
    plt.show()

print("detector locations: {}".format(Detectors))
print("Balloon location: {}".format(Balloon))
print("traveltimes: {}".format(traveltimes))
MiddleOfDetectors = np.array([0,0,0])
for Detector in Detectors:
    MiddleOfDetectors = MiddleOfDetectors + Detector/float(len(Detectors))

NumberOfDetectors = len(Detectors)
delta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))

#-------------------------------------------------------------------------------#
#                                   Fit n                                       #
#-------------------------------------------------------------------------------#
differences = np.zeros(len(indexofrefractionrange))

paths = []
times = []
distances = []


b_ballon = MiddleOfDetectors[2]
a_ballon = (Balloon[2]-b_ballon)/Balloon[0]

for number,n in enumerate(indexofrefractionrange):
    thetas = np.linspace(0,0.9,1000)
    delta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    correlation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    normedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    summedcorrelation = np.zeros(1000)

    for i in range(NumberOfDetectors):
        for j in range(NumberOfDetectors):
            if i < j:
                delta_t[i][j] = np.abs(traveltimes[i] - traveltimes[j])
                deltaz = np.linalg.norm(Detectors[i]-Detectors[j])
                delta_taccenten[i][j] = delta_taccent(thetas,np.abs(deltaz),n)
                correlation[i][j] = np.abs(delta_t[i][j] - delta_taccenten[i][j])
                normedcorrelation[i][j] = correlation[i][j]/np.trapz(correlation[i][j],thetas)

                summedcorrelation += normedcorrelation[i][j]


    angle_index = np.where(summedcorrelation == summedcorrelation.min())
    angle = thetas[angle_index] #zenith ofc

    a_planewave = np.tan(np.pi/2-angle)
    b_planewave = MiddleOfDetectors[2]
    angle_snell = np.arcsin(np.sin(angle)*n_icesurface)
    a_snell = np.tan(np.pi/2 - angle_snell)
    b_snell = -1*a_snell*(-1*b_planewave/a_planewave)
    XopBallonHoogte = (Balloon[2] - b_snell)/a_snell
    verschil = XopBallonHoogte - Balloon[0]

    differences[number] = np.abs(verschil)

#plot for the n that was found:
n_index = np.where(differences == differences.min())
n = indexofrefractionrange[n_index[0]]
if len(n) > 1:
    print("undetermined")
    print(n)
    n = n[0]
print(n)

thetas = np.linspace(0,0.9,1000)
delta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
correlation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
normedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
summedcorrelation = np.zeros(1000)

for i in range(NumberOfDetectors):
    for j in range(NumberOfDetectors):
        if i < j:
            delta_t[i][j] = np.abs(traveltimes[i] - traveltimes[j])
            deltaz = np.linalg.norm(Detectors[i]-Detectors[j])
            delta_taccenten[i][j] = delta_taccent(thetas,np.abs(deltaz),n)
            correlation[i][j] = np.abs(delta_t[i][j] - delta_taccenten[i][j])
            normedcorrelation[i][j] = correlation[i][j]/np.trapz(correlation[i][j],thetas)
            plt.plot(thetas,normedcorrelation[i][j])
            plt.show()

            summedcorrelation += normedcorrelation[i][j]

angle_index = np.where(summedcorrelation == summedcorrelation.min())
angle = thetas[angle_index] #zenith ofc
angle_snell = np.arcsin(np.sin(angle)*n_icesurface)

x = np.linspace(0,Balloon[0],1000)
a_refracted = np.tan(np.pi/2-angle)
b_refracted = MiddleOfDetectors[2]
a_snell = np.tan(np.pi/2-angle_snell)
b_snell = -1*a_snell*(-1*b_refracted/a_refracted)

y = a_refracted*x + b_refracted
indexwhensurface = np.where(y > 0)[0][0]
x = x[0:indexwhensurface]
y = y[0:indexwhensurface]
plt.plot(x,y,color="red",label="Plane Wave Reconstructed Path")
x_snell = np.linspace(x[-1],Balloon[0])
y = a_snell*x_snell + b_snell
plt.legend()
plt.plot(x_snell,y,color="red")

plt.show()

print("minimal difference between timing:")
print(differences.min())
if differences.min() > 1:
    print(f"{bcolors.FAIL}THIS FIT IS NOT USABLE{bcolors.ENDC}")
    sys.exit(1)
else:
    print(f"{bcolors.OKGREEN}This fit might be usable :) {bcolors.ENDC}")
n_index = np.where(differences == differences.min())
n_fit = indexofrefractionrange[n_index]
if len(n_fit) > 1:
    print("undetermined")
    print(n_fit)
    n_fit = n_fit[0]
print("index of refraction from fit at a depth of {}m : {}".format(MiddleOfDetectors[2],n_fit))
print("index of refraction from exponential model: {}".format(ice.get_index_of_refraction(MiddleOfDetectors)))
sys.exit(0)
