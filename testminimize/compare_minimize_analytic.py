from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import numpy as np
import logging
import configs


configmin = dict()
configmin['propagation'] = dict(
    attenuate_ice = True,
    focusing_limit = 2,
    focusing = False,
    radiopropa = dict(
        mode = 'minimizing',
        iter_steps_channel = [25., 2., .5], #unit is meter
        iter_steps_zenith = [.5, .05, .005], #unit is degree
        auto_step_size = False,
        max_traj_length = 10000) #unit is meter
)
configmin['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)


logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')
missed = 0
amountofvertices = 20

xpoints = np.random.uniform(low=-4000, high=4000.0, size=amountofvertices)
zpoints = np.random.uniform(low=-3000, high=0.0, size=amountofvertices)
points = [xpoints,zpoints]
final_point = np.array( [0., 0., -100.] ) * units.m

DirectDeltaTimes = []
RefractedDeltaTimes = [] 
ReflectedDeltaTimes = []

with alive_bar(amountofvertices,title='Calculating paths',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for i in range(amountofvertices):

        # Let's work in the y = 0 plane
        start_point = np.array([points[0][i],0.,points[1][i]])*units.m

        DirectDeltaTime = 0
        RefractedDeltaTime = 0 
        ReflectedDeltaTime = 0

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
            if SolType == 2:
                AmountOfRefracted += 1
                RefractedDeltaTime = prop.get_travel_time(Sol)
            if SolType == 3:
                AmountOfReflected += 1
                ReflectedDeltaTime = prop.get_travel_time(Sol)



        #analytic one:
        
        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configmin)

        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumberHybrid = prop.get_number_of_solutions()
        if (SolNumberHybrid == SolNumberAnalytic):
            for Sol in range(SolNumberHybrid):
                SolType = prop.get_solution_type(Sol) #1:direct,2:refracted and 3: reflected

                if SolType == 1:
                    AmountOfDirect -= 1
                    DirectDeltaTime -= prop.get_travel_time(Sol)
                if SolType == 2:
                    AmountOfRefracted -= 1
                    RefractedDeltaTime -= prop.get_travel_time(Sol)
                if SolType == 3:
                    AmountOfReflected -= 1
                    ReflectedDeltaTime -= prop.get_travel_time(Sol)


        if (AmountOfDirect == 0) and (AmountOfRefracted == 0) and (AmountOfReflected == 0):
            DirectDeltaTimes.append(DirectDeltaTime)
            RefractedDeltaTimes.append(RefractedDeltaTime)
            ReflectedDeltaTimes.append(ReflectedDeltaTime)
        else:
            missed += 1
            if (SolNumberHybrid == SolNumberAnalytic):
                print('prob \"non-correct identification\" on' + str(start_point))
            else:
                print('maybe a problem on' + str(start_point))

                
        bar()
print("there were {} \"special\" cases".format(missed))
plt.subplot(3, 1, 1)
plt.hist(DirectDeltaTimes/len(DirectDeltaTimes),bins=10)

plt.subplot(3, 1, 2)
plt.hist(RefractedDeltaTimes/len(RefractedDeltaTimes),bins=10)

plt.subplot(3, 1, 3)
plt.hist(ReflectedDeltaTimes/len(ReflectedDeltaTimes),bins=10)

plt.show()
