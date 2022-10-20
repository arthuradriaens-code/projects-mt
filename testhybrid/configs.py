from NuRadioMC.utilities import medium
from NuRadioReco.utilities import units

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
confighybrid['propagation']['radiopropa']['mode'] = 'hybrid minimizing'

