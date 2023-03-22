import numpy as np
import matplotlib.pyplot as plt
from NuRadioReco.modules.io.rno_g import rnogDataReader
from NuRadioReco.utilities import units
from NuRadioReco.modules import channelAddCableDelay
import logging
from NuRadioReco.detector import detector
from datetime import datetime
import NuRadioReco.modules.io.rno_g.readRNOGData
import radiotools.helper

#file location
filepath = "/home/arthur/Documents/thesis/data/interesting/station12/run687/combined.root"

#detector file path
pathToJson = "/home/arthur/Documents/thesis/programs/analysis-tools/rnog_analysis_tools/detector_json/RNO_season_2022.json"

#station number (id)
station_id = 12

AddCableDelay = channelAddCableDelay.channelAddCableDelay()
det = detector.Detector(json_filename = pathToJson)
det.update(datetime.now())

reader = NuRadioReco.modules.io.rno_g.readRNOGData.readRNOGData()
reader.begin(filepath)
n_channels = len([0,1,2,3])
time_differences = np.zeros((n_channels, n_channels))
reader.run(channels=np.array([0,1,2,3]),event_numbers={687:[6213]})
events = reader.get_events()
for event in events:
    print(event)
