from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig()

c = 299792458 #(m/s)
n = 1.78
v = c/n

logger = logging.getLogger('ray_tracing_modules')

# Let us work on the y = 0 plane
Detectors = np.zeros((4,3))
Detectors[0] = np.array([0., 0., -97.]) * units.m
Detectors[1] = np.array([0., 0., -96.]) * units.m
Detectors[2] = np.array([0., 0., -95.]) * units.m
Detectors[3] = np.array([0., 0., -94.]) * units.m
plt.plot([0,0,0,0], [-97,-96,-95,-94], 'bo')
traveltimes = []
def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degreetorad(deg):
    return deg*np.pi/180

paths = []
times = []
distances = []

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

plt.plot([33], [150], 'bo')

prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configh)
for detector in Detectors:
    start_point = np.array([ 33.0,0.,150.0])*units.m
    final_point = detector
    prop.set_start_and_end_point(start_point, final_point)
    prop.find_solutions()
    SolNumber = prop.get_number_of_solutions()
    for Sol in range(SolNumber):
        paths.append(prop.get_path(Sol))
        times.append(prop.get_travel_time(Sol))
        SolType = prop.get_solution_type(Sol)
        if SolType == 1:
            print('1 direct ray')
        if SolType == 2:
            print('1 refracted ray')
        if SolType == 3:
            print('1 reflected ray')
        print(CalcAngleToGround(prop.get_receive_vector(Sol)))
        x = np.linspace(paths[-1][0,0],start_point[0],1000)
        xlen = x[-1]-x[0]
        z = np.linspace(paths[-1][0,2],start_point[2],1000)
        zlen = z[-1]-z[0]
        diagonallen = np.sqrt(xlen*xlen + zlen*zlen) #(m)
        plt.plot(x,z,color='orange')
        traveltime = times[-1]/units.ns + (diagonallen/c)*(10**9) #ns
        traveltimes.append(traveltime)
        plt.plot(np.abs(paths[-1][:,0]),paths[-1][:,2],label="travel time = {0:.2f} nanoseconds".format(traveltime) ,color="orange")
plt.ylabel("vertical distance (m)")
plt.xlabel("horizontal distance (m)")
plt.title("Greenland simple trajectory with GL1 attenuation\n solved with hybrid ray tracer")
plt.show()

def delta_taccent(theta,deltaz):
    return ((np.cos(theta)*deltaz)/v)*(10**9)

thetas = np.linspace(0,np.pi,1000)

delta_taccenten01 = delta_taccent(thetas,1)
delta_t01 = traveltimes[0]-traveltimes[1]
plt.plot(thetas,np.abs(delta_t01 - delta_taccenten01))

delta_taccenten02 = delta_taccent(thetas,2)
delta_t02 = traveltimes[0]-traveltimes[2]
plt.plot(thetas,np.abs(delta_t02 - delta_taccenten02))

delta_taccenten03 = delta_taccent(thetas,3)
delta_t03 = traveltimes[0]-traveltimes[3]
plt.plot(thetas,np.abs(delta_t03 - delta_taccenten03))

delta_taccenten12 = delta_taccent(thetas,1)
delta_t12 = traveltimes[1]-traveltimes[2]
plt.plot(thetas,np.abs(delta_t12 - delta_taccenten12))

delta_taccenten13 = delta_taccent(thetas,2)
delta_t13 = traveltimes[1]-traveltimes[3]
plt.plot(thetas,np.abs(delta_t13 - delta_taccenten13))

delta_taccenten23 = delta_taccent(thetas,1)
delta_t23 = traveltimes[2]-traveltimes[3]
plt.plot(thetas,np.abs(delta_t23 - delta_taccenten23))

plt.xlabel("theta (rad)")
plt.ylabel("correlation")
plt.title("correlation versus zenith")

plt.show()
