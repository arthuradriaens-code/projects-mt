from __future__ import absolute_import, division, print_function
import numpy as np
from radiotools import helper as hp
from radiotools import plthelpers as php
from matplotlib import pyplot as plt
from NuRadioReco.utilities import units
from NuRadioMC.utilities import medium
from NuRadioMC.utilities import plotting
from six import iteritems
import h5py
import argparse
import json
import time
import os
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot NuRadioMC event list output.')
parser.add_argument('inputdir', type=str,
                    help='path to NuRadioMC hdf5 simulation outputs')
parser.add_argument('--trigger_name', type=str, default=None, nargs='+',
                    help='the name of the trigger that should be used for the plots')
parser.add_argument('--Veff', type=str,
                    help='specify json file where effective volume is saved as a function of energy')
args = parser.parse_args()

dirname = args.inputdir
plot_folder = 'plots'
if(not os.path.exists(plot_folder)):
    os.makedirs(plot_folder)

ReferenceArray = h5py.File(os.path.join(args.inputdir,'Exponential.hdf5'),'r')

for filename in os.listdir(args.inputdir):
    print(filename)
    if filename == "Exponential.hdf5":
        continue
    else:
        fin = h5py.File(os.path.join(args.inputdir,filename), 'r')
        print('the following triggeres where simulated: {}'.format(fin.attrs['trigger_names']))
        if(args.trigger_name is None):
            triggered = np.array(fin['triggered'])
            print("you selected any trigger")
            trigger_name = 'all'
        else:
            if(len(args.trigger_name) > 1):
                print("trigger {} selected which is a combination of {}".format(args.trigger_name[0], args.trigger_name[1:]))
                trigger_name = args.trigger_name[0]
                plot_folder = os.path.join(dirname, 'plots', filename, args.trigger_name[0])
                if(not os.path.exists(plot_folder)):
                    os.makedirs(plot_folder)
                triggered = np.zeros(len(fin['multiple_triggers'][:, 0]), dtype=np.bool)
                for trigger in args.trigger_name[1:]:
                    iTrigger = np.squeeze(np.argwhere(fin.attrs['trigger_names'] == trigger))
                    triggered = triggered | np.array(fin['multiple_triggers'][:, iTrigger], dtype=np.bool)
            else:
                trigger_name = args.trigger_name[0]
                iTrigger = np.argwhere(fin.attrs['trigger_names'] == trigger_name)
                triggered = np.array(fin['multiple_triggers'][:, iTrigger], dtype=np.bool)
                print("\tyou selected '{}'".format(trigger_name))
                plot_folder = os.path.join(dirname, 'plots', filename, trigger_name)
                if(not os.path.exists(plot_folder)):
                    os.makedirs(plot_folder)

        weights = np.array(fin['weights'])[triggered]
        n_events = fin.attrs['n_events']
        
        ###########################
        # plot neutrino direction
        ###########################

        """
        fig, ax = php.get_histogram((np.array(fin['zeniths'])[triggered])/ units.deg, weights=weights,
                                    ylabel='weighted entries', xlabel='zenith angle [deg]',
                                    bins=np.arange(0, 181, 5), figsize=(6, 6),title=str(filename))
        reffig, ax = php.get_histogram((np.array(ReferenceArray['zeniths'])[triggered] )/ units.deg, weights=weights,
                                    ylabel='weighted entries', xlabel='zenith angle [deg]',
                                    bins=np.arange(0, 181, 5), figsize=(6, 6),title=str(filename))
        diff=plt.bar(np.arange(0,181,45), 
             height=(fig-reffig), edgecolor='black', 
             linewidth=1.2, color='red',width = 1, align = 'edge') 

        plg.savefig(os.path.join(plot_folder, 'neutrino_direction_comparison_{}.png'.format(filename)))
        """
        fig, (ax1, ax2) = plt.subplots(1, 2)
        plot1=ax1.hist(np.array(fin['zeniths'])[triggered]/ units.deg, bins=np.arange(0,181,5),weights=weights,edgecolor='k',linewidth=1.0,color='blue')
        refplot=ax2.hist(np.array(ReferenceArray['zeniths'])[triggered] / units.deg, bins=np.arange(0,181,5),weights=weights,edgecolor='k',linewidth=1.0,color='blue')
        barxvalues = np.arange(0,181,5)
        barxvalues = np.delete(barxvalues, len(barxvalues)-1) #delete last element
        plt.clf()
        diff=plt.bar(barxvalues, 
                     height=(plot1[0]-refplot[0]), edgecolor='black', 
                     linewidth=1.2, color='red',width = 1, align = 'edge') 
        plt.title("{} - Exponential".format(filename[:-5]))
        plt.xlabel("zenith angle [deg]")
        plt.ylabel("weighted entries")
        plt.savefig(os.path.join(plot_folder,"{} - Exponential.pdf".format(filename[:-5])))


        czen = np.cos(np.array(fin['zeniths'])[triggered])
        ReferenceCzen = np.cos(np.array(ReferenceArray['zeniths'])[triggered])
        bins = np.linspace(-1, 1, 21)
        fig, ax = php.get_histogram(czen - ReferenceCzen, weights=weights,
                                    ylabel='weighted entries', xlabel='cos(zenith angle)',
                                    bins=bins, figsize=(6, 6),title=filename)

        # ax.set_xticks(np.arange(0, 181, 45))

        ax.set_title(trigger_name)
        fig.savefig(os.path.join(plot_folder, 'neutrino_direction_cos_comparison.png'))
