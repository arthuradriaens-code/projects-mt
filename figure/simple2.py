from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
import logging
logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')

# Let us work on the y = 0 plane
final_point4 = np.array( [0., 0., -95.] ) * units.m
rectangle4 = plt.Rectangle((-8,-100), 16, 10, fc='blue')
plt.gca().add_patch(rectangle4)

x1, y1 = [0, 0], [0, -95]
plt.plot(x1, y1, marker = '.', color='black')
plt.plot(800, -400, marker = 'x', color='red')
plt.axhline(y=0, color='grey',linestyle="--")

final_points = [final_point4]
def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degreetorad(deg):
    return deg*np.pi/180

AngleToFind = degreetorad(56)
print(AngleToFind)
paths = []
times = []
distances = []
start_point = np.array([-800,0.,-400.])*units.m
colors = ['orange','red','blue','green']

for i,final_point in enumerate(final_points): 
    prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
    prop.set_start_and_end_point(start_point, final_point)
    prop.find_solutions()
    SolNumber = prop.get_number_of_solutions()
    for Sol in range(SolNumber):
        paths.append(prop.get_path(Sol))
        times.append(prop.get_travel_time(Sol))
        plt.plot(np.abs(paths[-1][:,0]),paths[-1][:,2],label="travel time = {0:.2f} nanoseconds".format(times[-1]/units.ns),color=colors[i])


plt.ylabel("vertical distance (m)")
plt.xlabel("horizontal distance (m)")

#yticks = [-7,-20,-40,-60,-80,-100]
#plt.yticks(yticks)
#xticks = [0,20,40,60,80,100]
#plt.xticks(xticks)


plt.gca().invert_xaxis()
#plt.title("Greenland simple trajectory with GL1 attenuation\n solved with analytic ray tracer")
plt.savefig("detector-illustration.png", transparent=True)
plt.show()
