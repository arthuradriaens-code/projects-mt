from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
import os
from coordinate_system import CoordinateSystem
import datetime
import folium
from dateutil.tz import tzutc
import copy
import numpy as np
import csv


#gpx_file_name = 'SMT_20220630_112143.gpx'  #Jess
#gpx_file_name = 'SMT_20220701_231934.gpx' #Bob
#gpx_file_name = 'SMT_20220715_231621.gpx' #Bob
directory = os.fsencode('/home/arthur/Documents/thesis/data/sonde/gpx/')
with open('EventsBelow15Deg.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['Balloon filename (gpx)','station','expected timeframe'])
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        print(filename)
        if (filename[7] == "1") or (filename[7] == "3") or ((filename[9] != "9") and (filename[9] != "8") and (filename[9] != "7") and (filename[9] != "6")):
            print("not doing this one")
            continue
        gpx_file = os.fsdecode(directory) + filename
        i = 0
        track = 0
        for track in read_gpx_file(gpx_file):
            continue

        segment = track['segments'][0]

        stations = {21:[72.5874063909459,-38.4660301212611], 12:[72.6000868058195,-38.4962265332872],11:[72.5892267215905,-38.5022988244688], 13:[72.6109470001738,-38.4901465440588], 22:[72.598265271346,-38.4599355034766], 23:[72.6091242603966,-38.4538331609837] ,24:[72.6199833575357,-38.4477230792255]}
        coor = CoordinateSystem()
        for StationNumber in stations:
            StationCoordinates = stations[StationNumber]
            locallocationstation = coor.geodetic_to_enu(StationCoordinates[0],StationCoordinates[1],3251.9489147560234)
            for i in range(len(segment['lat'])):
                BalloonPosition = coor.geodetic_to_enu(segment['lat'][i],segment['lon'][i],segment['elevation'][i])
                HorizontalDistance = np.sqrt((BalloonPosition[0] - locallocationstation[0])**2 + (BalloonPosition[1] - locallocationstation[1])**2) 
                HeightDifference = BalloonPosition[2] + 95.5
                Angle = np.arctan(HorizontalDistance/HeightDifference)*180/np.pi

            #let's plot this
            the_map = create_folium_map(tiles='opentopomap')
            masked_segment = {key:[] for key in segment.keys()}
            masked_segment['elevation-up'] = segment['elevation-up']
            masked_segment['elevation-down'] = segment['elevation-down']

            for i, time in enumerate(segment['time']):
                for key in list(segment.keys())[:-2]:
                    masked_segment[key].append(segment[key][i])

            for key in list(segment.keys())[:-2]:
                if key in ['time','latlon']: continue
                masked_segment[key] = np.array(masked_segment[key])

            StartAndStopPositions = []
            StartAndStopTimes = []
            CloseEnough = False
            for i in range(len(segment['lat'])):
                BalloonPosition = coor.geodetic_to_enu(segment['lat'][i],segment['lon'][i],segment['elevation'][i])
                HorizontalDistance = np.sqrt((BalloonPosition[0] - locallocationstation[0])**2 + (BalloonPosition[1] - locallocationstation[1])**2) 
                HeightDifference = BalloonPosition[2] + 95.5
                Angle = np.arctan(HorizontalDistance/HeightDifference)*180/np.pi
                if abs(Angle) < 15:
                    StartAndStopPositions.append([segment['lat'][i],segment['lon'][i],segment['elevation'][i]])
                    StartAndStopTimes.append(segment['time'][i])
                    CloseEnough = True
                    
            if CloseEnough:
                print("The balloon from file {} gets close enough to detector {} from {} until {}".format(filename,StationNumber,StartAndStopTimes[0],StartAndStopTimes[-1]))

                add_segment_to_map(the_map, masked_segment,color_by='elevation')
    # Add, markers to the gps-locations of the detectors

                folium.Marker(
                    location=StationCoordinates,
                    popup="station {}".format(StationNumber),
                    icon=folium.Icon(color="blue"),
                ).add_to(the_map)

                folium.Marker(
                    location=[StartAndStopPositions[0][0],StartAndStopPositions[0][1]],
                    popup="Start of <15° at {}".format(StartAndStopTimes[0]),
                    icon=folium.Icon(color="green"),
                ).add_to(the_map)

                folium.Marker(
                    location=[StartAndStopPositions[-1][0],StartAndStopPositions[-1][1]],
                    popup="Stop of <15° at {}".format(StartAndStopTimes[-1]),
                    icon=folium.Icon(color="red"),
                ).add_to(the_map)

                boundary = the_map.get_bounds()
                the_map.fit_bounds(boundary, padding=(3, 3))
                the_map
    # To save the map:
                the_map.save('maps15/map_{}_{}.html'.format(filename,StationNumber))
                spamwriter.writerow([filename,StationNumber,'{} till {}'.format(StartAndStopTimes[0],StartAndStopTimes[-1])])
