from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import scipy.stats as st
import sys
import time
import numpy as np
import logging
import configs

# Save the current stdout so that we can revert sys.stdout after we complete
# our redirection
stdout_fileno = sys.stdout
StepsZenith = 0.7
SphereSize= 45

confighybrid = dict()
confighybrid['propagation'] = dict(
    attenuate_ice = True,
    focusing_limit = 2,
    focusing = False,
    radiopropa = dict(
        mode = 'hybrid minimizing',
        iter_steps_channel = [SphereSize, 2., .5], #unit is meter
        iter_steps_zenith = [StepsZenith, .05, .005], #unit is degree
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
        iter_steps_channel = [SphereSize, 2., .5], #unit is meter
        iter_steps_zenith = [StepsZenith, .05, .005], #unit is degree
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

xpoints = np.random.uniform(low=100, high=4000.0, size=amountofvertices)
zpoints = np.random.uniform(low=-3000, high=-100.0, size=amountofvertices) 
points = [xpoints,zpoints]
final_point = np.array( [0., 0., -100.] ) * units.m

DirectDeltaTimes = []
RefractedDeltaTimes = [] 
ReflectedDeltaTimes = []

DirectDeltaAZs = []
RefractedDeltaAZs = [] 
ReflectedDeltaAZs = []

GoodPointsItx = []
GoodPointsItz = []
BadPointsItx = []
BadPointsItz = []

GoodPointsHybx = []
GoodPointsHybz = []
BadPointsHybx = []
BadPointsHybz = []

TimeAn = []
TimeIt = []
TimeHyb = []


testcase = 0

with alive_bar(amountofvertices,title='Calculating paths using iterative raytracer',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for i in range(amountofvertices):

        # Let's work in the y = 0 plane
        start_point = np.array([points[0][i],0.,points[1][i]])*units.m

        DirectDeltaTime = 0
        RefractedDeltaTime = 0 
        ReflectedDeltaTime = 0

        DirectDeltaAZ = 0
        RefractedDeltaAZ = 0 
        ReflectedDeltaAZ = 0


        #iterative one:
        AmountOfDirect = 0
        AmountOfRefracted = 0
        AmountOfReflected = 0

        SecondOfSameRefracted = False
        
        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configit)

        prop.set_start_and_end_point(start_point, final_point)
        start = time.time()
        prop.find_solutions()
        end = time.time()
        TimeIt.append(end - start)
        SolNumberHybrid = prop.get_number_of_solutions()
        for Sol in range(SolNumberHybrid):
            SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

            if SolType == 1:
                AmountOfDirect += 1
                DirectDeltaTime = prop.get_travel_time(Sol)
                DirectDeltaAZ = ZenithAngle(prop.get_receive_vector(Sol))
            if SolType == 2:
                AmountOfRefracted += 1
                if AmountOfRefracted == 2:
                    RefractedDeltaTimeSec = prop.get_travel_time(Sol)
                    RefractedDeltaAZSec = ZenithAngle(prop.get_receive_vector(Sol))
                    SecondOfSameRefracted = True

                    #needed to sort values and compare right ones
                    RefractedDeltaTimeListIt = []
                    RefractedDeltaAZListIt = []
                    RefractedDeltaTimeListAn = []
                    RefractedDeltaAZListAn = []

                    RefractedDeltaTimeListIt.append(RefractedDeltaTime)
                    RefractedDeltaTimeListIt.append(RefractedDeltaTimeSec)
                    RefractedDeltaTimeArrayIt = np.array(RefractedDeltaTimeListIt)

                    RefractedDeltaAZListIt.append(RefractedDeltaAZ)
                    RefractedDeltaAZListIt.append(RefractedDeltaAZSec)
                    RefractedDeltaAZArrayIt = np.array(RefractedDeltaTimeListIt)
                else:
                    RefractedDeltaTime = prop.get_travel_time(Sol)
                    RefractedDeltaAZ = ZenithAngle(prop.get_receive_vector(Sol))
            if SolType == 3:
                AmountOfReflected += 1
                ReflectedDeltaTime = prop.get_travel_time(Sol)
                ReflectedDeltaAZ = ZenithAngle(prop.get_receive_vector(Sol))

        #analytic one:

        prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
        prop.set_start_and_end_point(start_point, final_point)
        start = time.time()
        prop.find_solutions()
        end = time.time()
        SolNumberAnalytic = prop.get_number_of_solutions()
        if (SolNumberHybrid == SolNumberAnalytic):
            for Sol in range(SolNumberAnalytic):
                SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

                if SolType == 1:
                    AmountOfDirect -= 1
                    DirectDeltaTime -= prop.get_travel_time(Sol)
                    if (np.abs(DirectDeltaTime/units.ns) > 1):
                        BadPointsItx.append(start_point[0])
                        BadPointsItz.append(start_point[2])
                        print(start_point)
                        print("is a bad point")
                        testcase += 1
                    else:
                        GoodPointsItx.append(start_point[0])
                        GoodPointsItz.append(start_point[2])
                    DirectDeltaAZ -= ZenithAngle(prop.get_receive_vector(Sol))
                    #DirectDeltaAZ = 0
                if SolType == 2:
                    AmountOfRefracted -= 1
                    if SecondOfSameRefracted:
                        RefractedDeltaTimeSec -= prop.get_travel_time(Sol)
                        RefractedDeltaAZSec -= ZenithAngle(prop.get_receive_vector(Sol))

                        #needed to sort values and compare right ones
                        RefractedDeltaTimeListAn.append(RefractedDeltaTimeSec)
                        RefractedDeltaTimeArrayAn = np.array(RefractedDeltaTimeListIt)

                        RefractedDeltaAZListAn.append(RefractedDeltaAZSec)
                        RefractedDeltaAZArrayAn = np.array(RefractedDeltaTimeListIt)
                    else:
                        RefractedDeltaTime -= prop.get_travel_time(Sol)
                        RefractedDeltaAZ -= ZenithAngle(prop.get_receive_vector(Sol))
                if SolType == 3:
                    AmountOfReflected -= 1
                    ReflectedDeltaTime -= prop.get_travel_time(Sol)
                    ReflectedDeltaAZ -= ZenithAngle(prop.get_receive_vector(Sol))
                    #ReflectedDeltaAZ = 0
        TimeAn.append((end - start))


        if (AmountOfDirect == 0) and (AmountOfRefracted == 0) and (AmountOfReflected == 0):

            DirectDeltaTimes.append(DirectDeltaTime/units.ns)
            ReflectedDeltaTimes.append(ReflectedDeltaTime/units.ns)
            DirectDeltaAZs.append(DirectDeltaAZ*(180/np.pi)*1000)
            ReflectedDeltaAZs.append(ReflectedDeltaAZ*(180/np.pi)*1000)
            if not SecondOfSameRefracted:
                RefractedDeltaTimes.append(RefractedDeltaTime/units.ns)
                RefractedDeltaAZs.append(RefractedDeltaAZ*(180/np.pi)*1000)
            else: 
                RefractedDeltaTimeArray = np.sort(RefractedDeltaTimeArrayIt/units.ns) - np.sort(RefractedDeltaTimeArrayAn/units.ns)
                RefractedDeltaAZArray = np.sort(RefractedDeltaAZArrayIt) - np.sort(RefractedDeltaAZArrayAn)
                for DeltaTime in RefractedDeltaTimeArray:
                    RefractedDeltaTimes.append(DeltaTime)
                for DeltaAZ in RefractedDeltaTimeArray:
                    RefractedDeltaAZs.append(DeltaAZ*(180/np.pi)*1000)
                    if np.abs(DeltaAZ*(180/np.pi)*1000) > 5000:
                        print("wtf on"+str(start_point))

                print("lesgoo")
        else:
            missed += 1
            if (SolNumberHybrid == SolNumberAnalytic):
                print('prob \"non-correct identification\" on' + str(start_point))
            else:
                print('maybe a problem on' + str(start_point))
                testcase += 1

                
        bar()
print("there were {} \"special\" cases".format(missed))
iterativemissed=testcase

print("the iterative raytracer time:" + str(np.sum(np.array(TimeIt))/amountofvertices) + "time an:" + str(np.sum(np.array(TimeAn))/amountofvertices) + "milliseconds")


plt.title("iterative")
plt.subplot(3, 2, 5)
hist,bin_edges = np.histogram(DirectDeltaTimes, bins=10000, range=None,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(DirectDeltaTimes,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="blue")
plt.ylabel("Direct rays")
plt.xlabel("nanoseconds")
plt.legend()

plt.subplot(3, 2, 3)
hist,bin_edges = np.histogram(RefractedDeltaTimes, bins=10000, range=None,weights=np.ones(len(RefractedDeltaTimes))/len(RefractedDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(RefractedDeltaTimes,weights=np.ones(len(RefractedDeltaTimes))/len(RefractedDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="orange")
plt.ylabel("Refracted rays")
plt.legend()

plt.subplot(3, 2, 1)
hist,bin_edges = np.histogram(ReflectedDeltaTimes, bins=10000,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(ReflectedDeltaTimes,weights=np.ones(len(ReflectedDeltaTimes))/len(ReflectedDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="green")
plt.ylabel("Reflected rays")
plt.title("arrival time difference")
plt.legend()

plt.subplot(3, 2, 2)
hist,bin_edges = np.histogram(DirectDeltaAZs, bins=10000, range=None,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(DirectDeltaAZs,weights=np.ones(len(DirectDeltaAZs))/len(DirectDeltaAZs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="green")
plt.title("arrival zenith difference")
plt.legend()

plt.subplot(3, 2, 4)
hist,bin_edges = np.histogram(RefractedDeltaAZs, bins=10000, range=None,weights=np.ones(len(RefractedDeltaAZs))/len(RefractedDeltaAZs), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(RefractedDeltaAZs,weights=np.ones(len(RefractedDeltaAZs))/len(RefractedDeltaAZs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="orange")
plt.legend()

plt.subplot(3, 2, 6)
hist,bin_edges = np.histogram(ReflectedDeltaAZs, bins=10000, range=None,weights=np.ones(len(ReflectedDeltaAZs))/len(ReflectedDeltaAZs), density=None)
percentile = Calc68(hist,bin_edges)
plt.hist(ReflectedDeltaAZs,weights=np.ones(len(ReflectedDeltaAZs))/len(ReflectedDeltaAZs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="blue")
plt.xlabel("millidegree")
plt.legend()

plt.savefig('iterative_comparison_N_{}.pdf'.format(amountofvertices),transparent=True)

#hybrid

DirectDeltaTimes = []
RefractedDeltaTimes = [] 
ReflectedDeltaTimes = []

DirectDeltaAZs = []
RefractedDeltaAZs = [] 
ReflectedDeltaAZs = []

GoodPointsx = []
GoodPointsz = []
BadPointsx = []
BadPointsz = []

with alive_bar(amountofvertices,title='Calculating paths using hybrid raytracer',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for i in range(amountofvertices):

        # Let's work in the y = 0 plane
        start_point = np.array([points[0][i],0.,points[1][i]])*units.m

        DirectDeltaTime = 0
        RefractedDeltaTime = 0 
        ReflectedDeltaTime = 0

        DirectDeltaAZ = 0
        RefractedDeltaAZ = 0 
        ReflectedDeltaAZ = 0


        #hybrid one:

        AmountOfDirect = 0
        AmountOfRefracted = 0
        AmountOfReflected = 0

        SecondOfSameRefracted = False
        
        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=confighybrid)

        prop.set_start_and_end_point(start_point, final_point)
        start = time.time()
        prop.find_solutions()
        end = time.time()
        SolNumberHybrid = prop.get_number_of_solutions()
        for Sol in range(SolNumberHybrid):
            SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

            if SolType == 1:
                AmountOfDirect += 1
                DirectDeltaTime = prop.get_travel_time(Sol)
                DirectDeltaAZ = ZenithAngle(prop.get_receive_vector(Sol))
            if SolType == 2:
                AmountOfRefracted += 1
                if AmountOfRefracted == 2:
                    RefractedDeltaTimeSec = prop.get_travel_time(Sol)
                    RefractedDeltaAZSec = ZenithAngle(prop.get_receive_vector(Sol))
                    SecondOfSameRefracted = True

                    #needed to sort values and compare right ones
                    RefractedDeltaTimeListIt = []
                    RefractedDeltaAZListIt = []
                    RefractedDeltaTimeListAn = []
                    RefractedDeltaAZListAn = []

                    RefractedDeltaTimeListIt.append(RefractedDeltaTime)
                    RefractedDeltaTimeListIt.append(RefractedDeltaTimeSec)
                    RefractedDeltaTimeArrayIt = np.array(RefractedDeltaTimeListIt)

                    RefractedDeltaAZListIt.append(RefractedDeltaAZ)
                    RefractedDeltaAZListIt.append(RefractedDeltaAZSec)
                    RefractedDeltaAZArrayIt = np.array(RefractedDeltaTimeListIt)
                else:
                    RefractedDeltaTime = prop.get_travel_time(Sol)
                    RefractedDeltaAZ = ZenithAngle(prop.get_receive_vector(Sol))
            if SolType == 3:
                AmountOfReflected += 1
                ReflectedDeltaTime = prop.get_travel_time(Sol)
                ReflectedDeltaAZ = ZenithAngle(prop.get_receive_vector(Sol))

        TimeHyb.append(end - start)
        
        #analytic one:

        prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumberAnalytic = prop.get_number_of_solutions()
        if (SolNumberHybrid == SolNumberAnalytic):
            for Sol in range(SolNumberAnalytic):
                SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

                if SolType == 1:
                    AmountOfDirect -= 1
                    DirectDeltaTime -= prop.get_travel_time(Sol)
                    if (np.abs(DirectDeltaTime/units.ns) > 1):
                        BadPointsHybx.append(start_point[0])
                        BadPointsHybz.append(start_point[2])
                        print(start_point)
                        print("is a bad point")
                    else:
                        GoodPointsHybx.append(start_point[0])
                        GoodPointsHybz.append(start_point[2])
                    DirectDeltaAZ -= ZenithAngle(prop.get_receive_vector(Sol))
                    #DirectDeltaAZ = 0
                if SolType == 2:
                    AmountOfRefracted -= 1
                    if SecondOfSameRefracted:
                        RefractedDeltaTimeSec -= prop.get_travel_time(Sol)
                        RefractedDeltaAZSec -= ZenithAngle(prop.get_receive_vector(Sol))

                        #needed to sort values and compare right ones
                        RefractedDeltaTimeListAn.append(RefractedDeltaTimeSec)
                        RefractedDeltaTimeArrayAn = np.array(RefractedDeltaTimeListIt)

                        RefractedDeltaAZListAn.append(RefractedDeltaAZSec)
                        RefractedDeltaAZArrayAn = np.array(RefractedDeltaTimeListIt)
                    else:
                        RefractedDeltaTime -= prop.get_travel_time(Sol)
                        RefractedDeltaAZ -= ZenithAngle(prop.get_receive_vector(Sol))
                if SolType == 3:
                    AmountOfReflected -= 1
                    ReflectedDeltaTime -= prop.get_travel_time(Sol)
                    ReflectedDeltaAZ -= ZenithAngle(prop.get_receive_vector(Sol))
                    #ReflectedDeltaAZ = 0



        if (AmountOfDirect == 0) and (AmountOfRefracted == 0) and (AmountOfReflected == 0):

            DirectDeltaTimes.append(DirectDeltaTime/units.ns)
            ReflectedDeltaTimes.append(ReflectedDeltaTime/units.ns)
            DirectDeltaAZs.append(DirectDeltaAZ*(180/np.pi)*1000)
            ReflectedDeltaAZs.append(ReflectedDeltaAZ*(180/np.pi)*1000)
            if not SecondOfSameRefracted:
                RefractedDeltaTimes.append(RefractedDeltaTime/units.ns)
                RefractedDeltaAZs.append(RefractedDeltaAZ*(180/np.pi)*1000)
            else: 
                RefractedDeltaTimeArray = np.sort(RefractedDeltaTimeArrayIt/units.ns) - np.sort(RefractedDeltaTimeArrayAn/units.ns)
                RefractedDeltaAZArray = np.sort(RefractedDeltaAZArrayIt) - np.sort(RefractedDeltaAZArrayAn)
                for DeltaTime in RefractedDeltaTimeArray:
                    RefractedDeltaTimes.append(DeltaTime)
                for DeltaAZ in RefractedDeltaTimeArray:
                    RefractedDeltaAZs.append(DeltaAZ*(180/np.pi)*1000)
                    if np.abs(DeltaAZ*(180/np.pi)*1000) > 5000:
                        print("wtf on"+str(start_point))
                print("lesgoo")
        else:
            missed += 1
            if (SolNumberHybrid == SolNumberAnalytic):
                print('prob \"non-correct identification\" on' + str(start_point))
            else:
                print('maybe a problem on' + str(start_point))

                
        bar()
print("there were {} \"special\" cases".format(missed))
hybridmissed=missed


print("the hybrid raytracer time:" + str(np.sum(np.array(TimeHyb))/amountofvertices) + "time an:" + str(np.sum(np.array(TimeAn))/amountofvertices) + "seconds")
sys.stdout = open('hybrid_sigmas', 'w')

plt.clf()
plt.title("hybrid")
plt.subplot(3, 2, 5)
hist,bin_edges = np.histogram(DirectDeltaTimes, bins=10000, range=None,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
sys.stdout.write(str(percentile))
sys.stdout.write("  ")
plt.hist(DirectDeltaTimes,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="blue")
plt.ylabel("Direct rays")
plt.xlabel("nanoseconds")
plt.legend()

plt.subplot(3, 2, 3)
hist,bin_edges = np.histogram(RefractedDeltaTimes, bins=10000, range=None,weights=np.ones(len(RefractedDeltaTimes))/len(RefractedDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
sys.stdout.write(str(percentile))
sys.stdout.write("  ")
plt.hist(RefractedDeltaTimes,weights=np.ones(len(RefractedDeltaTimes))/len(RefractedDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="orange")
plt.ylabel("Refracted rays")
plt.legend()

plt.subplot(3, 2, 1)
hist,bin_edges = np.histogram(ReflectedDeltaTimes, bins=10000,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
sys.stdout.write(str(percentile))
sys.stdout.write("  ")
plt.hist(ReflectedDeltaTimes,weights=np.ones(len(ReflectedDeltaTimes))/len(ReflectedDeltaTimes),bins=20,label="68%={0:.2f} ns".format(percentile),color="green")
plt.ylabel("Reflected rays")
plt.title("arrival time difference")
plt.legend()

plt.subplot(3, 2, 2)
hist,bin_edges = np.histogram(DirectDeltaAZs, bins=10000, range=None,weights=np.ones(len(DirectDeltaTimes))/len(DirectDeltaTimes), density=None)
percentile = Calc68(hist,bin_edges)
sys.stdout.write(str(percentile))
sys.stdout.write("  ")
plt.hist(DirectDeltaAZs,weights=np.ones(len(DirectDeltaAZs))/len(DirectDeltaAZs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="green")
plt.title("arrival zenith difference")
plt.legend()

plt.subplot(3, 2, 4)
hist,bin_edges = np.histogram(RefractedDeltaAZs, bins=10000, range=None,weights=np.ones(len(RefractedDeltaAZs))/len(RefractedDeltaAZs), density=None)
percentile = Calc68(hist,bin_edges)
sys.stdout.write(str(percentile))
sys.stdout.write("  ")
plt.hist(RefractedDeltaAZs,weights=np.ones(len(RefractedDeltaAZs))/len(RefractedDeltaAZs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="orange")
plt.legend()

plt.subplot(3, 2, 6)
hist,bin_edges = np.histogram(ReflectedDeltaAZs, bins=10000, range=None,weights=np.ones(len(ReflectedDeltaAZs))/len(ReflectedDeltaAZs), density=None)
percentile = Calc68(hist,bin_edges)
sys.stdout.write(str(percentile))
sys.stdout.write("  ")
plt.hist(ReflectedDeltaAZs,weights=np.ones(len(ReflectedDeltaAZs))/len(ReflectedDeltaAZs),bins=20,label="68%={0:.2f} mdeg".format(percentile),color="blue")
plt.xlabel("millidegree")
plt.legend()

plt.savefig('hybrid_comparison_N_{}.pdf'.format(amountofvertices),transparent=True)
sys.stdout.write(str(np.sum(np.array(TimeHyb))/np.sum(np.array(TimeAn))))
sys.stdout.write("  ")
sys.stdout.write(str(np.sum(np.array(TimeIt))/np.sum(np.array(TimeAn))))
sys.stdout.write("  ")
sys.stdout.write(str(np.sum(np.array(TimeHyb))/np.sum(np.array(TimeIt))))
sys.stdout.write("  ")

if (iterativemissed >= hybridmissed):
    sys.stdout.write("+")
else: 
    sys.stdout.write("X")

sys.stdout.write("\n")

# Close the file
sys.stdout.close()
# Restore sys.stdout to our old saved file handler
sys.stdout = stdout_fileno
