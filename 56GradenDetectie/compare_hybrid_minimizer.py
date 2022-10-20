from NuRadioMC.SignalProp import radioproparaytracing
from NuRadioMC.SignalProp import analyticraytracing
from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units
import matplotlib.pyplot as plt
import numpy as np
import logging
import time
from alive_progress import alive_bar


logging.basicConfig()

logger = logging.getLogger('ray_tracing_modules')
missed = 0

# Let us work on the y = 0 plane
final_point = np.array( [0., 0., -100.] ) * units.m

amountofvertices = 5
xpoints = np.random.uniform(low=-4000, high=4000.0, size=amountofvertices)
zpoints = np.random.uniform(low=-3000, high=0.0, size=amountofvertices)
points = [xpoints,zpoints]


hybridtimes = []
hybriddeltatimes = []
minimizetimes = []
minimizedeltatimes = []
iterativetimes = []
iterativedeltatimes = []
AnalyticTimes = []

#configs
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


configminim = dict()
configminim['propagation'] = dict(
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
configminim['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)
configminim['propagation']['radiopropa']['mode'] = 'minimizing'
 
confighyb = dict()
confighyb['propagation'] = dict(
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
confighyb['speedup'] = dict(
    delta_C_cut = 40 * units.degree
)

with alive_bar(amountofvertices,title='Calculating paths',length=20,bar='filling',spinner='dots_waves2') as bar: #fun progress bar, of course not needed for the program to function
    for i in range(amountofvertices):
        start_point = np.array([points[0][i],0.,points[1][i]])*units.m
        print('start point:'+str(start_point))
        start = time.time()

        paths = []
        times = []

        prop = analyticraytracing.ray_tracing(medium.greenland_simple(), attenuation_model='GL1')
        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumber = prop.get_number_of_solutions()
        for Sol in range(SolNumber):
            paths.append(prop.get_path(Sol))
            times.append(prop.get_travel_time(Sol))
        times.sort() 
        basepathtimes = times
         
        end = time.time()
        basetime = end - start
        AnalyticTimes.append(basetime)

        start = time.time()

        paths = []
        times = []

        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configminim)

        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumber = prop.get_number_of_solutions()
        for Sol in range(SolNumber):
            paths.append(prop.get_path(Sol))
            times.append(prop.get_travel_time(Sol))

        times.sort()
        if len(times) == len(basepathtimes):
            for i,pathtime in enumerate(times):
                minimizedeltatimes.append(pathtime-basepathtimes[i])
        else:
            missed += 1

        end = time.time()
        minimizetime = end - start
        minimizetimes.append(minimizetime)

        start = time.time()

        paths = []
        times = []

        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=confighyb)

        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumber = prop.get_number_of_solutions()
        for Sol in range(SolNumber):
            paths.append(prop.get_path(Sol))
            times.append(prop.get_travel_time(Sol))


        times.sort()
        if len(times) == len(basepathtimes):
            for i,pathtime in enumerate(times):
                hybriddeltatimes.append(pathtime-basepathtimes[i])
        else:
            missed += 1

        end = time.time()
        hybridtime = end - start
        hybridtimes.append(hybridtime)


        start = time.time()

        paths = []
        times = []

        prop = radioproparaytracing.radiopropa_ray_tracing(medium.greenland_simple(), attenuation_model='GL1',config=configit)

        prop.set_start_and_end_point(start_point, final_point)
        prop.find_solutions()
        SolNumber = prop.get_number_of_solutions()
        for Sol in range(SolNumber):
            paths.append(prop.get_path(Sol))
            times.append(prop.get_travel_time(Sol))

        times.sort()
        if len(times) == len(basepathtimes):
            for i,pathtime in enumerate(times):
                iterativedeltatimes.append(pathtime-basepathtimes[i])
        else:
            missed += 1

        end = time.time()
        iterativetime = end - start
        iterativetimes.append(iterativetime)
        bar()

totaltimes = amountofvertices 
deltatotaltimes = totaltimes - missed
print("average hybrid time is {} and nanosecond error {}".format(np.sum(np.array(hybridtimes))/totaltimes,np.sum(np.array(hybriddeltatimes))/deltatotaltimes))
print("average minimize time is {} and nanosecond error {}".format(np.sum(np.array(minimizetimes))/totaltimes, np.sum(np.array(minimizedeltatimes))/deltatotaltimes))
print("average iterative time is {} and nanosecond error {}".format(np.sum(np.array(iterativetimes))/totaltimes,np.sum(np.array(iterativedeltatimes))/deltatotaltimes))
print("average analytic time is {}".format(np.sum(np.array(AnalyticTimes ))/totaltimes))
