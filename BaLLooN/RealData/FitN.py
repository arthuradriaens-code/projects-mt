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
import scipy.optimize as opt
from scipy.signal import find_peaks
from scipy.stats import norm
from astropy import modeling

station_id = 21
event_id = 117
channel_ids = [5,7]
n_channels = len(channel_ids)
prefix="/mnt/usb"
gpx_file = prefix+"/sonde/gpx/SMT_20220726_111605.gpx"
GivenTime = "2022/07/26/11/18/41"
rootfile = prefix+"/RNO-G-DATA/station21/run1441/combined.root"
n_cut = [1.665,1.675] #if 2 or more peaks are observed will cut n within this border

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

def S(t,sampling_rate,amplitude):
    A = amplitude
    f = 0.403
    return A*np.sin(2*np.pi*f*t)
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
def InSumErrorOnIndex(thetamin,deltaz,C):
    #timely error is assumed in ns and deltaz in meters
    c = 0.2997925 #(m/ns)
    return ((c/(deltaz*np.cos(thetamin)))**2)*(1+C**2)
def ErrorOnIndex(epsilon,delta_t,summederror):
    return 2*(1+epsilon*0.01)*delta_t*np.sqrt(summederror)
def SimpleError(epsilon,delta_t,Delta_z,summedcorr,thetamin):
    c = 0.2997925 #(m/ns)
    return 2*c*(1+epsilon*0.01)*delta_t/Delta_z*np.sqrt(summedcorr**2+1)/np.cos(thetamin)
def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]
def FWHM(X,Y):
    hmx = half_max_x(X,Y)
    FWHM = hmx[1] - hmx[0]
    return FWHM
    
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
indexofrefractionrange = np.linspace(1.45,1.9,100000)
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
FoundEventId=False
for event in data_reader.run():
    if event.get_id()==event_id:
        FoundEventId=True
        break
if not FoundEventId:
    print(f"{bcolors.FAIL}EVENT ID NOT FOUND{bcolors.ENDC}")
    sys.exit(1)

# get detector information at the specified event time
det = NuRadioReco.detector.detector.Detector(json_filename="/home/arthur/Universiteit/master-proef/analysis-tools/rnog_analysis_tools/detector_json/RNO_season_2022.json", 
                                             antenna_by_depth=False)
station = event.get_station(station_id)
det.update(station.get_station_time())
passband = np.array([0.15,0.60])*units.GHz


locallocationstation = np.array(det.get_absolute_position(station_id))

# Ok now we have our balloon position and the detector position, 
# We'll set the BalloonPosition relative to the detector:

Balloon = np.array([BalloonPosition[0]-locallocationstation[0],BalloonPosition[1] - locallocationstation[1],BalloonPosition[2]] - locallocationstation[2])

#And let's rotate to only have an x component, no y:
r = np.sqrt(Balloon[0]**2 + Balloon[1]**2)
Balloon[1] = 0
Balloon[0] = r


#-------------------------------------------------------------------------------#
#                           Get correlation with sine                           #
#-------------------------------------------------------------------------------#

Detectors = []
corrs = []
cable_delays = np.zeros(n_channels)

for i,channel_id in enumerate(channel_ids):
    channel = station.get_channel(channel_id)
    target_sampling_rate = 10.0*units.GHz
    channel.resample(target_sampling_rate)
    channel = station.get_channel(channel_id)
    cable_delays[i] = det.get_cable_delay(station_id,channel_id)
    Detectors.append(np.array(det.get_relative_position(station_id,channel_id)) * units.m)
    channel_voltages = channel.get_filtered_trace(passband, 'butter', 6) 
    channel_sample_rate = channel.get_sampling_rate()/units.GHz
    channel_spectrum = channel.get_frequency_spectrum()
    channel_frequencies = channel.get_frequencies()
    #plt.axvline(x=0.403,linestyle='dashed',color="grey")
    #plt.xlabel("Frequency (GHz)")
    #plt.ylabel("Counts")
    #plt.title("Frequency spectrum of channel {}".format(channel_id))
    #plt.plot(channel_frequencies,channel_spectrum)
    #plt.show()
    channel_times = channel.get_times()
    #plt.plot(channel_times,channel_voltages,label="measured data")
    lmin,lmax = hl_envelopes_idx(channel_voltages)
    #plt.xlabel("time (nanoseconds)")
    #plt.ylabel("Voltage")
    #plt.title("Voltage i.f.o time for channel {}".format(channel_id))
    #plt.legend()
    #plt.show()
    t = np.arange(0,3/(0.403*units.GHz),1/channel_sample_rate)
    amplitude = 0.007
    template = S(t,channel_sample_rate,amplitude)
    #get correlation sine and signal
    corr = radiotools.helper.get_normalized_xcorr(
        channel_voltages,
        template
    )
    corrs.append(corr)
    times = np.linspace(0,len(corr)*(1/channel_sample_rate),len(corr))
    #plt.plot(times,corr)
    #plt.xlabel("time (nanoseconds)")
    #plt.ylabel("correlation")
    #plt.title("correlation i.f.o delta time for channel {}".format(channel_id))
    #plt.show()

#-------------------------------------------------------------------------------#
#                           Get expected time differences                       #
#-------------------------------------------------------------------------------#

# This is needed to limit the range in which we'll look

simtimes = []
simtraveltimes = []
simpaths = []

configh = dict()
configh['propagation'] = dict(
    attenuate_ice = True,
    focusing_limit = 2,
    focusing = False,
    radiopropa = dict(
        mode = 'hybrid minimizing',
        iter_steps_channel = [25., 2., .5], #unit is meter
        iter_steps_zenith = [.5, .05, .005], #unit is degree
        auto_step_size = False,
        max_traj_length = 10000) #unit is meter
)
configh['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)
prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configh)
for detector in Detectors:
    start_point = Balloon
    final_point = detector
    prop.set_start_and_end_point(start_point, final_point)
    prop.find_solutions()
    SolNumber = prop.get_number_of_solutions()
    for Sol in range(SolNumber):
        simpaths.append(prop.get_path(Sol))
        simtimes.append(prop.get_travel_time(Sol))
        SolType = prop.get_solution_type(Sol)
        x = np.linspace(simpaths[-1][0,0],start_point[0],1000)
        xlen = x[-1]-x[0]
        z = np.linspace(simpaths[-1][0,2],start_point[2],1000)
        zlen = z[-1]-z[0]
        diagonallen = np.sqrt(xlen*xlen + zlen*zlen) #(m)
        simtraveltime = simtimes[-1]/units.ns + (diagonallen/c)*(10**9) #ns
        simtraveltimes.append(simtraveltime)

expected_delta_t = np.zeros((n_channels, n_channels))
for i in range(n_channels):
    for j in range(i+1,n_channels):
        expected_delta_t[i][j] = simtraveltimes[i]-simtraveltimes[j]
        expected_delta_t[j][i] = simtraveltimes[j]-simtraveltimes[i]

#-------------------------------------------------------------------------------#
#                           Get observed time differences                       #
#-------------------------------------------------------------------------------#
# NOTE: This doesn't seem to work yet

   
#correlate the correlations, maximum will be offset
delta_t = np.zeros((n_channels, n_channels))
FoundPeak = False
for i in range(n_channels):
    for j in range(i+1,n_channels):
        corr = radiotools.helper.get_normalized_xcorr(
                corrs[i],
                corrs[j] 
            )
        t_offsets = (np.arange(
                    -len(corr) // 2,
                    len(corr) // 2, dtype=int
                ) / target_sampling_rate) - (cable_delays[i] - cable_delays[j])
        peaks, _ = find_peaks(corr, height=0)
        #plt.plot(t_offsets,corr)
        #plt.plot(t_offsets[peaks],corr[peaks],"x")
        #plt.xlabel("time (nanoseconds)")
        #plt.ylabel("correlation")
        #plt.title("correlation of sine correlations of channels {} and {}".format(channel_ids[i],channel_ids[j]))
        #plt.show()
        peaktimes = t_offsets[peaks]
        for peaktime in peaktimes:
            if np.abs(peaktime - expected_delta_t[i][j]) < 1.24: 
                FoundPeak=True
                print("difference between channels {} and {} is {}ns".format(channel_ids[i],channel_ids[j],peaktime))
                print("expected difference between channels {} and {} is {}ns".format(channel_ids[i],channel_ids[j],expected_delta_t[i][j]))
                delta_t[i, j] = peaktime
                delta_t[j, i] = -peaktime
        if not FoundPeak:
            print(f"{bcolors.FAIL}EITHER THE SIGNAL FROM CHANNEL {channel_ids[i]} OR {channel_ids[j]} ISN'T USABLE{bcolors.ENDC}")
            sys.exit(1)


        #delta_t[i, j] = t_offsets[np.argmax(corr)]
        #delta_t[j, i] = -t_offsets[np.argmax(corr)]

MiddleOfDetectors = np.array([0,0,0])
for Detector in Detectors:
    MiddleOfDetectors = MiddleOfDetectors + Detector/float(len(Detectors))

NumberOfDetectors = len(Detectors)
#-------------------------------------------------------------------------------#
#                               Find epsilon                                    #
#-------------------------------------------------------------------------------#

def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degrees(SomethingInRads):
    return SomethingInRads*180/np.pi
def delta_taccent(theta,deltaz,n):
    v = c/n
    return ((np.cos(theta)*deltaz)/v)*(10**9)


Epsindexofrefractionrange = np.linspace(1.6,1.9,5000)

Epstraveltimes = []
Epspaths = []
Epstimes = []
Epsdistances = []

ice = medium.greenland_simple()
prop = radioproparaytracing.radiopropa_ray_tracing(ice, attenuation_model='GL1',config=configh)
for detector in Detectors:
    start_point = Balloon
    final_point = detector
    prop.set_start_and_end_point(start_point, final_point)
    prop.find_solutions()
    SolNumber = prop.get_number_of_solutions()
    for Sol in range(SolNumber):
        Epspaths.append(prop.get_path(Sol))
        Epstimes.append(prop.get_travel_time(Sol))
        x = np.linspace(Epspaths[-1][0,0],start_point[0],1000)
        xlen = x[-1]-x[0]
        z = np.linspace(Epspaths[-1][0,2],start_point[2],1000)
        zlen = z[-1]-z[0]
        diagonallen = np.sqrt(xlen*xlen + zlen*zlen) #(m)
        Epstraveltime = Epstimes[-1]/units.ns + (diagonallen/c)*(10**9) #ns
        Epstraveltimes.append(Epstraveltime)


Epsdifferences = np.zeros(len(Epsindexofrefractionrange))

for number,n in enumerate(Epsindexofrefractionrange):
    #find plane wave
    Epsthetas = np.linspace(0,0.9,1000)
    NumberOfDetectors = len(Detectors)
    Epsdelta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))
    Epsdelta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    Epscorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    Epsnormedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
    Epssummedcorrelation = np.zeros(1000)

    for i in range(NumberOfDetectors):
        for j in range(NumberOfDetectors):
            if i < j:
                Epsdelta_t[i][j] = np.abs(Epstraveltimes[i] - Epstraveltimes[j])
                Epsdeltaz = np.linalg.norm(Detectors[i]-Detectors[j])
                Epsposition = (Detectors[0] + Detectors[1])/2
                Epsdelta_taccenten[i][j] = delta_taccent(Epsthetas,np.abs(Epsdeltaz),n)
                Epscorrelation[i][j] = np.abs(Epsdelta_t[i][j] - Epsdelta_taccenten[i][j])
                Epsnormedcorrelation[i][j] = Epscorrelation[i][j]/np.trapz(Epscorrelation[i][j],Epsthetas)

                Epssummedcorrelation += Epsnormedcorrelation[i][j]

    Epsangle_index = np.where(Epssummedcorrelation == Epssummedcorrelation.min())
    Epsangle = Epsthetas[Epsangle_index] #zenith ofc
    b_ballon = MiddleOfDetectors[2]
    a_ballon = (Balloon[2]-b_ballon)/Balloon[0]
    a_planewave = np.tan(np.pi/2-Epsangle)
    b_planewave = MiddleOfDetectors[2]

    angle_snell = np.arcsin(np.sin(Epsangle)*1.27) 
    a_snell = np.tan(np.pi/2-angle_snell)
    b_snell = -1*a_snell*(-1*b_planewave/a_planewave)
    XopBallonHoogte = (Balloon[2] - b_snell)/a_snell
    verschil = XopBallonHoogte - Balloon[0]

    Epsdifferences[number] = np.abs(verschil)

Epsn_index = np.where(Epsdifferences == Epsdifferences.min())
Epsn_fit = Epsindexofrefractionrange[Epsn_index]
if len(Epsn_fit) > 1:
    Epsn_fit = Epsn_fit[0]
Epsdirect_angle = np.pi/2 - np.arctan(a_ballon)
n_actual = ice.get_index_of_refraction(MiddleOfDetectors)

Epsilon = (100*(Epsn_fit - n_actual)/n_actual)

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
correctedindices = indexofrefractionrange/(1+Epsilon*0.01)
inversediffs = 1/differences

plt.plot(correctedindices,inversediffs)
plt.xlabel("index of refraction")
plt.ylabel("1/distance to balloon  ($m^{-1}$)")
plt.show()

if n_cut:
    #quick and dirty, could be faster using numpy
    for index,n in enumerate(correctedindices):
        if n > n_cut[0]:
            indexleft = index
            break
    for index,n in enumerate(correctedindices[::-1]):
        if n < n_cut[1]:
            indexright = len(correctedindices) - index
            break
    correctedindices = correctedindices[indexleft:indexright]
    inversediffs = inversediffs[indexleft:indexright]

plt.plot(correctedindices,inversediffs,label="data")
plt.xlabel("index of refraction")
plt.ylabel("1/distance to balloon  ($m^{-1}$)")
#fitter = modeling.fitting.LevMarLSQFitter()
#model = modeling.models.Gaussian1D()   # depending on the data you need to give some initial values
#fitted_model = fitter(model, correctedindices, inversediffs)
#plt.plot(correctedindices, fitted_model(correctedindices),label="gaussian fit")
peaks, _ = find_peaks(inversediffs, height=0)
n_FittedIndex = peaks[int(len(peaks)/2)]
n_fitted = correctedindices[n_FittedIndex]
plt.plot(n_fitted,inversediffs[n_FittedIndex],marker="o",label="fitted index of refraction")
print("fitted n index:")
print(n_FittedIndex)
for index in range(int(len(correctedindices)/2)):
    total = np.trapz(inversediffs,correctedindices)
    partial = np.trapz(inversediffs[n_FittedIndex-index:n_FittedIndex+index],correctedindices[n_FittedIndex-index:n_FittedIndex+index])
    if partial/total > 0.68:
        n_error68 = correctedindices[n_FittedIndex+index] - correctedindices[n_FittedIndex]
        n_bounds68 = [n_fitted - correctedindices[n_FittedIndex+index] + correctedindices[n_FittedIndex],n_fitted + correctedindices[n_FittedIndex+index] - correctedindices[n_FittedIndex]]
        break

for index in range(int(len(correctedindices)/2)):
    total = np.trapz(inversediffs,correctedindices)
    partial = np.trapz(inversediffs[n_FittedIndex-index:n_FittedIndex+index],correctedindices[n_FittedIndex-index:n_FittedIndex+index])
    if partial/total > 0.95:
        n_error95 = correctedindices[n_FittedIndex+index] - correctedindices[n_FittedIndex]
        n_bounds95 = [n_fitted - correctedindices[n_FittedIndex+index] + correctedindices[n_FittedIndex],n_fitted + correctedindices[n_FittedIndex+index] - correctedindices[n_FittedIndex]]
        break
plt.title("mean = {0:.4f}, 68% error: {1:.4f} and 95% error: {2:.4f}".format(n_fitted,n_error68,n_error95))
plt.vlines(n_bounds68,0,inversediffs[n_FittedIndex],colors="orange",label="68% bounds")
plt.vlines(n_bounds95,0,inversediffs[n_FittedIndex],colors="red",label="95% bounds")
plt.legend()
plt.show()

UnderSRError = 0

n = n_fitted

thetas = np.linspace(0,0.9,1000)
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

            #error -->
            angle_index_err = np.where(normedcorrelation[i][j] == normedcorrelation[i][j].min()) 
            angle_err = thetas[angle_index] #zenith ofc 

            UnderSRError += InSumErrorOnIndex(angle_err,deltaz,normedcorrelation[i][j].min())
            #<--

            #plt.plot(thetas,normedcorrelation[i][j])
            #plt.show()

            summedcorrelation += normedcorrelation[i][j]


angle_index = np.where(summedcorrelation == summedcorrelation.min())
angle = thetas[angle_index] #zenith ofc

Error = ErrorOnIndex(Epsilon,1/target_sampling_rate,UnderSRError)

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

print("minimal difference in meters between plane wave and actual:")
print(differences.min())
if differences.min() > 1:
    print(f"{bcolors.FAIL}THIS FIT IS NOT USABLE{bcolors.ENDC}")
    sys.exit(1)
else:
    print(f"{bcolors.OKGREEN}This fit might be usable :) {bcolors.ENDC}")
print("depth = {}m".format(MiddleOfDetectors[2]))
print("index of refraction from exponential model: {}".format(ice.get_index_of_refraction(MiddleOfDetectors)))
print("Epsilon: {}%".format(Epsilon))
print("index of refraction from fit after correction with epsilon: {} \u00B1 {}".format(n_fitted,n_error68))
sys.exit(0)
