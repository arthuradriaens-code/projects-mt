from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')

# Let us work on the y = 0 plane
final_point = np.array( [0., 0., -100.] ) * units.m

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

prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configh)

start_point = np.array([ 482.80938104,0.,-105.50417894])*units.m
prop.set_start_and_end_point(start_point, final_point)
prop.find_solutions()
SolNumber = prop.get_number_of_solutions()
for Sol in range(SolNumber):
    paths.append(prop.get_path(Sol))
    times.append(prop.get_travel_time(Sol))
    plt.plot(np.abs(paths[-1][:,0]),paths[-1][:,2],label="travel time = {0:.2f} nanoseconds".format(times[-1]/units.ns))
    SolType = prop.get_solution_type(Sol)
    if SolType == 1:
        print('1 direct ray')
    if SolType == 2:
        print('1 refracted ray')
    if SolType == 3:
        print('1 reflected ray')
    print(CalcAngleToGround(prop.get_receive_vector(Sol)))

plt.ylabel("vertical distance (m)")
plt.xlabel("horizontal distance (m)")

plt.gca().invert_xaxis()
plt.title("Greenland simple trajectory with GL1 attenuation\n solved with hybrid ray tracer")

plt.legend()
plt.show()
