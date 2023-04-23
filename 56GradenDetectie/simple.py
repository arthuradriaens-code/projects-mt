from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
import logging
import time
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
prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
start_point = np.array([800,0.,-405.50417894])*units.m
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
plt.ylabel("vertical distance (m)")
plt.xlabel("horizontal distance (m)")

xendnobend = 90
rico_a_NoBend = (100 - 7)/(xendnobend)
rico_b_NoBend = -100
x_No = np.arange(0,xendnobend,0.1)
y_No = rico_a_NoBend*x_No + rico_b_NoBend

plt.gca().invert_xaxis()
plt.title("Greenland simple trajectory with GL1 attenuation\n solved with analytic ray tracer")

plt.legend()
plt.show()
