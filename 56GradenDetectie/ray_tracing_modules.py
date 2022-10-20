from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')

# Let us work on the y = 0 plane
final_point = np.array( [0., 0., -100.] ) * units.m
xcoordinates = -np.arange(0,100,0.1)

def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degreetorad(deg):
    return deg*np.pi/180

anglestoground = [0,5,10,20,40]

for angletoground in anglestoground:
    AngleToFind = degreetorad(56-angletoground)
    print(AngleToFind)
    paths = []
    times = []
    distances = []

    for xcoordinate in xcoordinates:
        prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
        start_point = np.array([xcoordinate,0.,-7.])*units.m
        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumber = prop.get_number_of_solutions()
        for Sol in range(SolNumber):
            LaunchVector = prop.get_launch_vector(Sol)
            Angle = CalcAngleToGround(LaunchVector)
            if AngleToFind-0.01 < Angle < AngleToFind+0.01:
                print('found!')
                paths.append(prop.get_path(Sol))
                times.append(prop.get_travel_time(Sol))
                distances.append(-xcoordinate)
        prop = prop.reset_solutions()

    print("The travel time is:" + str(times[0]/units.s * units.ns))
    print("The distance to the origin is:" + str(distances[0]))
    plt.plot(np.abs(paths[-1][:,0]),paths[-1][:,2],label="zenith angle {0}°\ntravel time = {1:.2f} nanoseconds".format(angletoground,times[0]/units.ns))
    plt.ylabel("vertical distance (m)")
    plt.xlabel("horizontal distance (m)")

    xendnobend = 90/np.tan(degreetorad(34))
    rico_a_NoBend = (-100 +7)/(xendnobend)
    rico_b_NoBend = -7
    x_No = np.arange(0,xendnobend,0.1)
    y_No = rico_a_NoBend*x_No + rico_b_NoBend
    #plt.plot(x_No,y_No,color='grey',label="No bending",linestyle="--")

yticks = [-7,-20,-40,-60,-80,-100]
plt.yticks(yticks)
xticks = [0,20,40,60,80,100]
plt.xticks(xticks)
plt.gca().invert_xaxis()
plt.title("Greenland simple trajectory with GL1 attenuation")

plt.legend()
plt.show()
