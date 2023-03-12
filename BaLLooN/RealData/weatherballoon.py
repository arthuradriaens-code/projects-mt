from gpxplotter import read_gpx_file, create_folium_map, add_segment_to_map
import datetime
import folium
from dateutil.tz import tzutc
import copy
import numpy as np

the_map = create_folium_map(tiles='opentopomap')

#gpx_file_name = 'SMT_20220630_112143.gpx'  #Jess
#gpx_file_name = 'SMT_20220701_231934.gpx' #Bob
#gpx_file_name = 'SMT_20220715_231621.gpx' #Bob
gpx_file_name = '/media/RNO-G/data/sonde/gpx/SMT_20230107_231633.gpx'
i = 0
for track in read_gpx_file(gpx_file_name):
    continue

track.keys()

segment = track['segments'][0]
segment.keys()

segment['time'][0]

segment['time'][-1]

#launch_time = datetime.datetime(2022,6,30,11,5,0,0,tzutc()) #Jess
#launch_time = datetime.datetime(2022,7,1,23,19,0,0,tzutc()) #Bob
#launch_time

masked_segment = {key:[] for key in segment.keys()}
masked_segment['elevation-up'] = segment['elevation-up']
masked_segment['elevation-down'] = segment['elevation-down']
for i, time in enumerate(segment['time']):

    for key in list(segment.keys())[:-2]:
        masked_segment[key].append(segment[key][i])

for key in list(segment.keys())[:-2]:
    if key in ['time','latlon']: continue
    masked_segment[key] = np.array(masked_segment[key])

add_segment_to_map(the_map, masked_segment,color_by='elevation')
# Add, markers to the gps-locations of the detectors

folium.Marker(
    location=[72.5874063909459,-38.4660301212611],
    popup="Amaroq (station 21)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.598265271346,-38.4599355034766],
    popup="Avinngaq (station 22)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.6091242603966,-38.4538331609837],
    popup="Ukaliatsiaq (station 23)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.6199833575357,-38.4477230792255],
    popup="Qappik (station 24)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.5892267215905,-38.5022988244688],
    popup="Nanoq (station 11)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.6000868058195,-38.4962265332872],
    popup="Terianniaq (station 12)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.6109470001738,-38.4901465440588],
    popup="Ukaleq (station 13)",
    icon=folium.Icon(color="green"),
).add_to(the_map)

folium.Marker(
    location=[72.58279265212887, -38.45581495328228],
    popup="DISC",
    icon=folium.Icon(color="green"),
).add_to(the_map)


boundary = the_map.get_bounds()
the_map.fit_bounds(boundary, padding=(3, 3))
the_map
# To save the map:
the_map.save('map_000.html')
