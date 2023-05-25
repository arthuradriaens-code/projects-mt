from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import logging
logging.basicConfig()
from alive_progress import alive_bar

def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degrees(SomethingInRads):
    return SomethingInRads*180/np.pi
def delta_taccent(theta,deltaz,n):
    v = c/n
    return ((np.cos(theta)*deltaz)/v)*(10**9)
def eta(n,delta_t_observed,theta_b,delta_z):
    delta_t_computed = delta_taccent(theta_b,delta_z,n)
    return delta_t_observed - delta_t_computed
def f(n,NumberOfDetectors,theta_b,delta_t,delta_z):
    fvar = 0
    for i in range(NumberOfDetectors):
        for j in range(NumberOfDetectors):
            if i < j:
                fvar += eta(n,delta_t[i][j],theta_b,delta_z[i][j])**2
    return fvar
 

c = 299792458 #(m/s)
ice = medium.greenland_simple()

logger = logging.getLogger('ray_tracing_modules')


# Let us work on the y = 0 plane
Detectors = np.zeros((4,3))
Detectors[0] = np.array([0., 0., -97.]) * units.m
Detectors[1] = np.array([0., 0., -96.]) * units.m
Detectors[2] = np.array([0., 0., -95.]) * units.m
Detectors[3] = np.array([0., 0., -94.]) * units.m

configh = dict()
configh['propagation'] = dict(
    attenuate_ice = True,
    focusing_limit = 2,
    focusing = False,
    radiopropa = dict(
        mode = 'hybrid minimizing',
        iter_steps_channel = [25., 2., .5, .05,0.005], #unit is meter
        iter_steps_zenith = [.5, .05, .005, .001,0.0001], #unit is degree
        auto_step_size = False,
        max_traj_length = 10000) #unit is meter
)
configh['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)

zcoordinates = np.linspace(100,600,5)
xcoordinates = np.linspace(6,800,30)
indexofrefractionrange = np.linspace(1.27,2.5,5000)

with alive_bar(len(zcoordinates),title='Calculating relative angles',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for zcoordinate in zcoordinates:
        RelativeAccuracy = np.zeros(len(xcoordinates))
        BalloonAngle = np.zeros(len(xcoordinates))
        for s,xcoordinate in enumerate(xcoordinates):
            traveltimes = []
            paths = []
            times = []
            distances = []
            Balloon = np.array([xcoordinate,0.,500.0])*units.m
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

            NumberOfDetectors = len(Detectors)
            delta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))
            delta_z = np.zeros((NumberOfDetectors,NumberOfDetectors))

            for i in range(NumberOfDetectors):
                for j in range(NumberOfDetectors):
                    if i < j:
                        delta_t[i][j] = traveltimes[i] - traveltimes[j]
                        delta_z[i][j] = np.linalg.norm(Detectors[i]-Detectors[j])
            
            position = np.array([0,0,-95.5])
            n_actual = ice.get_index_of_refraction(position)
            print("actual n:")
            print(n_actual)
            b_ballon = -95.5
            a_ballon = (Balloon[2]-b_ballon)/Balloon[0]
            theta_b = np.pi/2 - np.arctan(a_ballon)
            
            root = optimize.minimize(f,x0=n_actual-0.1,args=(NumberOfDetectors,theta_b,delta_t,delta_z),method='Nelder-Mead')

            n_fit = root.x
            if len(n_fit) > 1:
                print("undetermined")
                print(n_fit)
                n_fit = n_fit[0]

            print("fitted n:")
            print(n_fit)

            BalloonAngle[s] = degrees(theta_b)
            RelativeAccuracy[s] = 100*(n_fit - n_actual)/n_actual

        bar()
            
        plt.plot(BalloonAngle,RelativeAccuracy,label="height = {}".format(zcoordinate))
        plt.xlabel("Direct angle (°)")
        plt.ylabel("$\epsilon$ (%)")

plt.legend()
plt.title("$\epsilon$ in function of the Direct angle")
plt.show()
