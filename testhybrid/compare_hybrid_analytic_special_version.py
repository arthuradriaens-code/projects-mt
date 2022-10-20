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
        mode = 'iterative',
        iter_steps_channel = [25., 2., .5], #unit is meter
        iter_steps_zenith = [.5, .05, .005], #unit is degree
        auto_step_size = False,
        max_traj_length = 10000) #unit is meter
)
confighybrid['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)

def ZenithAngle(a):
    lena = np.sqrt(np.dot(a,a)) #normalize
    return np.arccos(np.dot(a,np.array([0,0,1]))/lena)

logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')
missed = 0
amountofvertices = 40

xpoints = np.random.uniform(low=-4000, high=4000.0, size=amountofvertices)
zpoints = np.random.uniform(low=-3000, high=0.0, size=amountofvertices)
points = [xpoints,zpoints]
final_point = np.array( [0., 0., -100.] ) * units.m

DirectDeltaTimes = []
RefractedDeltaTimes = [] 
ReflectedDeltaTimes = []

DirectDeltaAZs = []
RefractedDeltaAZs = [] 
ReflectedDeltaAZs = []


with alive_bar(amountofvertices,title='Calculating paths',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for i in range(amountofvertices):

        # Let's work in the y = 0 plane
        start_point = np.array([points[0][i],0.,points[1][i]])*units.m

        DirectDeltaTimeAn = []
        RefractedDeltaTimeAn = [] 
        ReflectedDeltaTimeAn = []

        DirectDeltaAZAn = [] 
        RefractedDeltaAZAn = [] 
        ReflectedDeltaAZAn = []

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
                DirectDeltaTimeAn.append(prop.get_travel_time(Sol))
                DirectDeltaAZAn.append(ZenithAngle(prop.get_receive_vector(Sol)))
            if SolType == 2:
                AmountOfRefracted += 1
                RefractedDeltaTimeAn.append(prop.get_travel_time(Sol))
                RefractedDeltaAZAn.append(ZenithAngle(prop.get_receive_vector(Sol)))
            if SolType == 3:
                AmountOfReflected += 1
                ReflectedDeltaTimeAn.append(prop.get_travel_time(Sol))
                ReflectedDeltaAZAn.append(ZenithAngle(prop.get_receive_vector(Sol)))

        DirectDeltaTimeAn = np.sort(DirectDeltaTimeAn)
        DirectDeltaAZAn = np.sort(DirectDeltaAZAn)

        RefractedDeltaTimeAn = np.sort(RefractedDeltaTimeAn)
        RefractedDeltaAZAn = np.sort(RefractedDeltaAZAn)

        ReflectedDeltaTimeAn = np.sort(ReflectedDeltaTimeAn)
        ReflectedDeltaAZAn = np.sort(ReflectedDeltaAZAn)

        #hybrid one: 
        DirectDeltaTimeHyb = []
        RefractedDeltaTimeHyb = [] 
        ReflectedDeltaTimeHyb = []

        DirectDeltaAZHyb = [] 
        RefractedDeltaAZHyb = [] 
        ReflectedDeltaAZHyb = []

       
        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=confighybrid)

        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumberHybrid = prop.get_number_of_solutions()
        if (SolNumberHybrid == SolNumberAnalytic):
            for Sol in range(SolNumberHybrid):
                SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

                if SolType == 1:
                    AmountOfDirect -= 1
                    DirectDeltaTimeHyb.append(prop.get_travel_time(Sol))
                    DirectDeltaAZHyb.append(ZenithAngle(prop.get_receive_vector(Sol)))
                if SolType == 2:
                    AmountOfRefracted -= 1
                    RefractedDeltaTimeHyb.append(prop.get_travel_time(Sol))
                    RefractedDeltaAZHyb.append(ZenithAngle(prop.get_receive_vector(Sol)))
                if SolType == 3:
                    AmountOfReflected -= 1
                    ReflectedDeltaTimeHyb.append(prop.get_travel_time(Sol))
                    ReflectedDeltaAZHyb.append(ZenithAngle(prop.get_receive_vector(Sol)))


        DirectDeltaTimeHyb = np.sort(DirectDeltaTimeHyb)
        DirectDeltaAZHyb = np.sort(DirectDeltaAZHyb)

        RefractedDeltaTimeHyb = np.sort(RefractedDeltaTimeHyb)
        RefractedDeltaAZHyb = np.sort(RefractedDeltaAZHyb)

        ReflectedDeltaTimeHyb = np.sort(ReflectedDeltaTimeHyb)
        ReflectedDeltaAZHyb = np.sort(ReflectedDeltaAZHyb)



        if (AmountOfDirect == 0) and (AmountOfRefracted == 0) and (AmountOfReflected == 0):
            for i in range(len(DirectDeltaTimeHyb)):
                DirectDeltaTime = DirectDeltaTimeAn[i] - DirectDeltaTimeHyb[i] 
                DirectDeltaTimes.append(DirectDeltaTime/units.ns)
            for i in range(len(RefractedDeltaTimeHyb)):
                RefractedDeltaTime = RefractedDeltaTimeAn[i] - RefractedDeltaTimeHyb[i] 
                RefractedDeltaTimes.append(RefractedDeltaTime/units.ns)
            for i in range(len(ReflectedDeltaTimeHyb)):
                ReflectedDeltaTime = ReflectedDeltaTimeAn[i] - ReflectedDeltaTimeHyb[i] 
                ReflectedDeltaTimes.append(ReflectedDeltaTime/units.ns)
            for i in range(len(DirectDeltaAZHyb)):
                DirectDeltaAZ = DirectDeltaAZAn[i] - DirectDeltaAZHyb[i] 
                DirectDeltaAZs.append(DirectDeltaAZ*(180/np.pi)*1000)
            for i in range(len(RefractedDeltaAZHyb)):
                RefractedDeltaAZ = RefractedDeltaAZAn[i] - RefractedDeltaAZHyb[i] 
                RefractedDeltaAZs.append(RefractedDeltaAZ*(180/np.pi)*1000)
            for i in range(len(ReflectedDeltaAZHyb)):
                ReflectedDeltaAZ = ReflectedDeltaAZAn[i] - ReflectedDeltaAZHyb[i] 
                ReflectedDeltaAZs.append(ReflectedDeltaAZ*(180/np.pi)*1000)

        else:
            missed += 1
            if (SolNumberHybrid == SolNumberAnalytic):
                print('prob \"non-correct identification\" on' + str(start_point))
            else:
                print('maybe a problem on' + str(start_point))

                
        bar()
print("there were {} \"special\" cases".format(missed))

#plots
plt.subplot(3, 2, 5)
hist,bins = np.histogram(DirectDeltaTimes)
quantile1 = np.quantile(DirectDeltaTimes,.16)
quantile2 = np.quantile(DirectDeltaTimes,.84)
quantile = quantile2-quantile1
weights = np.ones_like(DirectDeltaTimes)/float(len(DirectDeltaTimes))
plt.hist(DirectDeltaTimes,weights=weights,bins=20,label="68%={0:.2f} ns".format(quantile),color="blue")
plt.ylabel("Direct rays")
plt.xlabel("nanoseconds")
plt.legend()

plt.subplot(3, 2, 3)
hist,bins = np.histogram(RefractedDeltaTimes)
quantile1 = np.quantile(RefractedDeltaTimes,.16)
quantile2 = np.quantile(RefractedDeltaTimes,.84)
quantile = quantile2-quantile1
weights = np.ones_like(RefractedDeltaTimes)/float(len(RefractedDeltaTimes))
plt.hist(RefractedDeltaTimes,weights=weights,bins=20,label="68%={0:.2f} ns".format(quantile),color="orange")
plt.ylabel("Refracted rays")
plt.legend()

plt.subplot(3, 2, 1)
hist,bins = np.histogram(ReflectedDeltaTimes)
quantile1 = np.quantile(ReflectedDeltaTimes,.16)
quantile2 = np.quantile(ReflectedDeltaTimes,.84)
quantile = quantile2-quantile1
weights = np.ones_like(ReflectedDeltaTimes)/float(len(ReflectedDeltaTimes))
rv = st.rv_discrete(values=(ReflectedDeltaTimes, weights/weights.sum()))
quantileinterval = rv.interval(0.68)
quantile = quantileinterval[1]-quantileinterval[0]
plt.hist(ReflectedDeltaTimes,weights=weights,bins=20,label="68%={0:.2f} ns".format(quantile),color="green")
plt.ylabel("Reflected rays")
plt.title("arrival time difference")
plt.legend()

plt.subplot(3, 2, 2)
hist,bins = np.histogram(ReflectedDeltaAZs)
quantile1 = np.quantile(ReflectedDeltaAZs,.16)
quantile2 = np.quantile(ReflectedDeltaAZs,.84)
quantile = quantile2-quantile1
weights = np.ones_like(ReflectedDeltaAZs)/float(len(ReflectedDeltaAZs))
plt.hist(ReflectedDeltaAZs,weights=weights,bins=20,label="68%={0:.2f} mdeg".format(quantile))
plt.title("arrival zenith difference")
plt.legend()

plt.subplot(3, 2, 4)
hist,bins = np.histogram(RefractedDeltaAZs)
quantile1 = np.quantile(RefractedDeltaAZs,.16)
quantile2 = np.quantile(RefractedDeltaAZs,.84)
quantile = quantile2-quantile1 
weights = np.ones_like(RefractedDeltaAZs)/float(len(RefractedDeltaAZs))
plt.hist(RefractedDeltaAZs,weights=weights,bins=20,label="68%={0:.2f} mdeg".format(quantile),color="orange")
plt.legend()

plt.subplot(3, 2, 6)
hist,bins = np.histogram(DirectDeltaAZs)
quantile1 = np.quantile(DirectDeltaAZs,.16)
quantile2 = np.quantile(DirectDeltaAZs,.84)
quantile = quantile2-quantile1 
weights = np.ones_like(DirectDeltaAZs)/float(len(DirectDeltaAZs))
plt.hist(DirectDeltaAZs,weights=weights,bins=20,label="68%={0:.2f} mdeg".format(quantile),color="blue")
plt.xlabel("millidegree")
plt.legend()

plt.savefig('iterative_comparison_N_{}.pdf'.format(amountofvertices),transparent=True)
plt.show()
