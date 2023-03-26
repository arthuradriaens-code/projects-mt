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
import time_difference_functions

#file location
filepath = "/mnt/usb/RNO-G-DATA/station23/run691/combined.root"

#detector file path
pathToJson = "/mnt/usb/detector_json/RNO_season_2022.json"

#station number (id)
station_id = 23
channel_1 = 5
channel_2 = 6
template = None
det = detector.Detector(json_filename=pathToJson)

test = time_difference_functions.get_time_differences(station_id,channel_1,channel_2,det,passband=None)
print(test)
