from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import scipy.stats as st
import numpy as np
import logging
import configs


confighybrid = dict()
confighybrid['propagation'] = dict(
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
confighybrid['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)

configit = dict()
configit['propagation'] = dict(
    attenuate_ice = True,
    focusing_limit = 2,
    focusing = False,
    radiopropa = dict(
        mode = 'iterative',
        iter_steps_channel = [25., 2., .5], #unit is meter
        iter_steps_zenith = [.5, .05, .005], #unit is degree
        auto_step_size = False,
        max_traj_length = 10000) #unit is meter
)
configit['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)


def ZenithAngle(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,1]))/lena)

def Calc68(hist,bins):
    totalleft = 0
    totalright = 0
    for i,histpoint in enumerate(hist):
        totalleft += histpoint
        if totalleft >= 0.16:
            leftindex = i
            break
    for i,histpoint in enumerate(hist[::-1]):
        totalright += histpoint
        if totalright >= 0.16:
            rightindex = i
            break
    rightindex = len(hist)-rightindex-1
    while hist[leftindex:rightindex].sum() >= 0.68:
        leftindex += 1
    leftindex -=1
    while hist[leftindex:rightindex].sum() >= 0.68:
        rightindex -= 1
    rightindex +=1
    return (bins[rightindex] - bins[leftindex])

logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')
missed = 0
amountofvertices = 1000

xpoints = np.random.uniform(low=-4000, high=4000.0, size=amountofvertices)
zpoints = np.random.uniform(low=-3000, high=0.0, size=amountofvertices)
points = [xpoints,zpoints]
final_point = np.array( [0., 0., -100.] ) * units.m

DirectDeltaTimes = []
RefractedDeltaTimes = [] 
ReflectedDeltaTimes = []

DirectDeltaVAs = []
RefractedDeltaVAs = [] 
ReflectedDeltaVAs = []

DirectDeltaAZs = []
RefractedDeltaAZs = [] 
ReflectedDeltaAZs = []


with alive_bar(amountofvertices,title='Calculating paths',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for i in range(amountofvertices):

        # Let's work in the y = 0 plane
        start_point = np.array([points[0][i],0.,points[1][i]])*units.m

        DirectDeltaTime = 0
        RefractedDeltaTime = 0 
        ReflectedDeltaTime = 0

        DirectDeltaVA = 0 
        RefractedDeltaVA = 0 
        ReflectedDeltaVA = 0


        DirectDeltaAZ = 0
        RefractedDeltaAZ = 0 
        ReflectedDeltaAZ = 0

        SecondOfSameRefracted = False

        #analytic one:

        AmountOfDirect = 0
        AmountOfRefracted = 0
        AmountOfReflected = 0

        prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumberAnalytic = prop.get_number_of_solutions()
        for Sol in range(SolNumberAnalytic):
            SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

            if SolType == 1:
                AmountOfDirect += 1
                DirectDeltaTime = prop.get_travel_time(Sol)
                DirectDeltaVA = ZenithAngle(prop.get_receive_vector(Sol))
                #DirectDeltaAZ = 0
            if SolType == 2:
                AmountOfRefracted += 1
                if AmountOfRefracted == 2:
                    RefractedDeltaTimeSec = prop.get_travel_time(Sol)
                    RefractedDeltaVASec = ZenithAngle(prop.get_receive_vector(Sol))
                    SecondOfSameRefracted = True
                else:
                    RefractedDeltaTime = prop.get_travel_time(Sol)
                    RefractedDeltaVA = ZenithAngle(prop.get_receive_vector(Sol))
                #RefractedDeltaAZ = 0 
            if SolType == 3:
                AmountOfReflected += 1
                ReflectedDeltaTime = prop.get_travel_time(Sol)
                ReflectedDeltaVA = ZenithAngle(prop.get_receive_vector(Sol))
                #ReflectedDeltaAZ = 0

        #hybrid one:
        
        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=confighybrid)

        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumberHybrid = prop.get_number_of_solutions()
        if (SolNumberHybrid == SolNumberAnalytic):
            for Sol in range(SolNumberHybrid):
                SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

                if SolType == 1:
                    AmountOfDirect -= 1
                    DirectDeltaTime -= prop.get_travel_time(Sol)
                    DirectDeltaVA -= ZenithAngle(prop.get_receive_vector(Sol))
                if SolType == 2:
                    AmountOfRefracted -= 1
                    if AmountOfRefracted == 1:
                        RefractedDeltaTimeSec -= prop.get_travel_time(Sol)
                        RefractedDeltaVASec -= ZenithAngle(prop.get_receive_vector(Sol))
                    else:
                        RefractedDeltaTime -= prop.get_travel_time(Sol)
                        RefractedDeltaVA -= ZenithAngle(prop.get_receive_vector(Sol))
                if SolType == 3:
                    AmountOfReflected -= 1
                    ReflectedDeltaTime -= prop.get_travel_time(Sol)
                    ReflectedDeltaVA -= ZenithAngle(prop.get_receive_vector(Sol))

        if (AmountOfDirect == 0) and (AmountOfRefracted == 0) and (AmountOfReflected == 0):
            DirectDeltaTimes.append(DirectDeltaTime/units.ns)
            RefractedDeltaTimes.append(RefractedDeltaTime/units.ns)
            ReflectedDeltaTimes.append(ReflectedDeltaTime/units.ns)

            DirectDeltaVAs.append(DirectDeltaVA*(180/np.pi)*1000)
            RefractedDeltaVAs.append(RefractedDeltaVA*(180/np.pi)*1000)
            ReflectedDeltaVAs.append(ReflectedDeltaVA*(180/np.pi)*1000)
            if SecondOfSameRefracted: 
                RefractedDeltaTimes.append(RefractedDeltaTimeSec/units.ns)
                RefractedDeltaVAs.append(RefractedDeltaVASec*(180/np.pi)*1000)
                print("lesgoo")
        else:
            missed += 1
            if (SolNumberHybrid == SolNumberAnalytic):
                print('prob \"non-correct identification\" on' + str(start_point))
            else:
                print('maybe a problem on' + str(start_point))

                
        bar()
print("there were {} \"special\" cases".format(missed))

plt.subplot(3, 2, 5)
hist,bin_edges = np.histogram(DirectDeltaTimes, bins=5000, range=None,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(DirectDeltaTimes,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="blue")
plt.ylabel("Direct rays")
plt.xlabel("nanoseconds")
plt.legend()

plt.subplot(3, 2, 3)
hist,bin_edges = np.histogram(RefractedDeltaTimes, bins=5000, range=None,weights=np.ones(len(RefractedDeltaTimes))/len(RefractedDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(RefractedDeltaTimes,weights=np.ones(len(RefractedDeltaTimes))/len(RefractedDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="orange")
plt.ylabel("Refracted rays")
plt.legend()

plt.subplot(3, 2, 1)
hist,bin_edges = np.histogram(ReflectedDeltaTimes, bins=5000,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(ReflectedDeltaTimes,weights=np.ones(len(ReflectedDeltaTimes))/len(ReflectedDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="green")
plt.ylabel("Reflected rays")
plt.title("arrival time difference")
plt.legend()

plt.subplot(3, 2, 2)
hist,bin_edges = np.histogram(DirectDeltaVAs, bins=5000, range=None,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(DirectDeltaVAs,weights=np.ones(len(DirectDeltaVAs))/len(DirectDeltaVAs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="green")
plt.title("arrival zenith difference")
plt.legend()

plt.subplot(3, 2, 4)
hist,bin_edges = np.histogram(RefractedDeltaVAs, bins=5000, range=None,weights=np.ones(len(RefractedDeltaVAs))/len(RefractedDeltaVAs), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(RefractedDeltaVAs,weights=np.ones(len(RefractedDeltaVAs))/len(RefractedDeltaVAs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="orange")
plt.legend()

plt.subplot(3, 2, 6)
hist,bin_edges = np.histogram(ReflectedDeltaVAs, bins=5000, range=None,weights=np.ones(len(ReflectedDeltaVAs))/len(ReflectedDeltaVAs), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(ReflectedDeltaVAs,weights=np.ones(len(ReflectedDeltaVAs))/len(ReflectedDeltaVAs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="blue")
plt.xlabel("millidegree")
plt.legend()

plt.savefig('hybrid_comparison_N_{}.pdf'.format(amountofvertices),transparent=True)
plt.show()
