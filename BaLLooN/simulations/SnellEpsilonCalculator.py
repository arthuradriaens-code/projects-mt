from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import numpy as np
import logging
logging.basicConfig()
from alive_progress import alive_bar

def CalcAngleToGround(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,-1]))/lena)
def degreetorad(deg):
    return deg*np.pi/180

c = 299792458 #(m/s)
ice = medium.greenland_simple()

logger = logging.getLogger('ray_tracing_modules')

def degrees(SomethingInRads):
    return SomethingInRads*180/np.pi

# Let us work on the y = 0 plane
Detectors = np.zeros((2,3))
Detectors[0] = np.array([0,0,-95]) * units.m
Detectors[1] = np.array([0.,0.,-96]) * units.m
#Detectors[2] = np.array([0., 0., -95.631]) * units.m
#Detectors[3] = np.array([0., 0., -94.615]) * units.m

MiddleOfDetectors = np.array([0,0,0])
for Detector in Detectors:
    MiddleOfDetectors = MiddleOfDetectors + Detector/float(len(Detectors))

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

xcoordinates = [96]
indexofrefractionrange = np.linspace(1.6,1.9,5000)
RelativeAccuracy = np.zeros(len(xcoordinates))
BalloonAngle = np.zeros(len(xcoordinates))

with alive_bar(len(xcoordinates),title='Calculating relative angles',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for s,xcoordinate in enumerate(xcoordinates):
        traveltimes = []
        paths = []
        times = []
        distances = []
        Balloon = np.array([xcoordinate,0.,1000])*units.m
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

        def delta_taccent(theta,deltaz,n):
            v = c/n
            return ((np.cos(theta)*deltaz)/v)*(10**9)

        differences = np.zeros(len(indexofrefractionrange))

        for number,n in enumerate(indexofrefractionrange):
            #find plane wave
            thetas = np.linspace(0,0.9,1000)
            NumberOfDetectors = len(Detectors)
            delta_t = np.zeros((NumberOfDetectors,NumberOfDetectors))
            delta_taccenten = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
            correlation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
            normedcorrelation = np.zeros((NumberOfDetectors,NumberOfDetectors,1000))
            summedcorrelation = np.zeros(1000)

            for i in range(NumberOfDetectors):
                for j in range(NumberOfDetectors):
                    if i < j:
                        delta_t[i][j] = np.abs(traveltimes[i] - traveltimes[j])
                        deltaz = np.linalg.norm(Detectors[i]-Detectors[j])
                        position = (Detectors[0] + Detectors[1])/2
                        delta_taccenten[i][j] = delta_taccent(thetas,np.abs(deltaz),n)
                        correlation[i][j] = np.abs(delta_t[i][j] - delta_taccenten[i][j])
                        normedcorrelation[i][j] = correlation[i][j]/np.trapz(correlation[i][j],thetas)

                        summedcorrelation += normedcorrelation[i][j]

            angle_index = np.where(summedcorrelation == summedcorrelation.min())
            angle = thetas[angle_index] #zenith ofc
            b_ballon = MiddleOfDetectors[2]
            a_ballon = (Balloon[2]-b_ballon)/Balloon[0]
            a_planewave = np.tan(np.pi/2-angle)
            b_planewave = MiddleOfDetectors[2]

            angle_snell = np.arcsin(np.sin(angle)*1.27) 
            a_snell = np.tan(np.pi/2-angle_snell)
            b_snell = -1*a_snell*(-1*b_planewave/a_planewave)
            XopBallonHoogte = (Balloon[2] - b_snell)/a_snell
            verschil = XopBallonHoogte - Balloon[0]

            differences[number] = np.abs(verschil)

        n_index = np.where(differences == differences.min())
        n_fit = indexofrefractionrange[n_index]
        if len(n_fit) > 1:
            print("undetermined")
            print(n_fit)
            n_fit = n_fit[0]
        print(n_fit)

        direct_angle = np.pi/2 - np.arctan(a_ballon)

        n_actual = ice.get_index_of_refraction(position)
        print(n_actual)

        BalloonAngle[s] = degrees(direct_angle)
        print("Epsilon: {}".format(100*(n_fit - n_actual)/n_actual))

        bar()
