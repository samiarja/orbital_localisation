from typing import Final
import httplib2
import requests
import os
import numpy as np
from tqdm import tqdm
from datetime import datetime, timezone
from bs4 import BeautifulSoup, SoupStrainer

'''
https://scihub.copernicus.eu/twiki/do/view/SciHubUserGuide/FullTextSearch?redirectedfrom=SciHubUserGuide.3FullTextSearch
https://scihub.copernicus.eu/dhus/#/home
'''

position_x      = []
position_y      = []
position_z      = []
earth_diameter = 6371
START  = '2022-01-24T20:07:00'
FINISH = '2022-01-24T20:15:00'
url = 'https://data.nasa.gov/browse?q=ISS+COORDS&sortBy=newest&utf8=%E2%9C%93&page='

if not os.path.isdir("ISS_TOPO"):
    os.mkdir("ISS_TOPO")

date_START = datetime.strptime(START, "%Y-%m-%dT%H:%M:%S")

def fetchData():
    match=[]
    with open('ISS_TOPO/ISS_COORDS_URL.txt') as file:
        lines = [line.rstrip() for line in file]
        for num,line in enumerate(lines):
            if ("https://data.nasa.gov/Space-Science/ISS_COORDS_{}-{}-").format("%02d" % date_START.year, "%02d" % date_START.month) in line:
                match += [line]
        for num,line in enumerate(match):
            day_value = line[55:-10]
            if day_value < "%02d" % date_START.day:
                match = line
                return line[47:-10]

def downloadData(link):
    r = requests.get('https://nasa-public-data.s3.amazonaws.com/iss-coords/' + link + '/ISS_OEM/ISS.OEM_J2K_EPH.txt', allow_redirects=True)
    open('ISS_TOPO/ISS.OEM_J2K_EPH_' + link + '.txt', 'wb').write(r.content)
    with open('ISS_TOPO/ISS.OEM_J2K_EPH_' + link + '.txt') as file:
        lines = [line.rstrip() for line in file]
    return lines

def findTimeRange(trajectory_data):
    with open(trajectory_data) as file:
        lines = [line.rstrip() for line in file]

def webcraping(url):
    f = open("ISS_TOPO/ISS_COORDS_URL.txt", "w")
    http = httplib2.Http()
    for page in tqdm(range(1, 16)):
        status, response = http.request((url + '{}').format(page))
        responseurl = requests.get((url + '{}').format(page))
        if responseurl.status_code == 200:
            for link in BeautifulSoup(response, parse_only=SoupStrainer('a')):
                if link.has_attr('href'):
                    data = link['href']
                    if 'https://data.nasa.gov/Space-Science/' in data:
                        f.write(link['href'])
                        f.write("\n")
        else:
            f.close()
    f.close()

def headerFinder(infile):
    for num,line in enumerate(infile):
        if line.startswith(START + '.000'):
            start_time = num
        if line.startswith(FINISH + '.000'):
            finish_time = num
    return start_time, finish_time

def interpolate(x1,x2):
    position_interpolation = np.linspace(float(x1), float(x2), num=4)
    return position_interpolation

def xyz_to_lat_lon(x,y,z):
    lat = np.degrees(np.arcsin(z/earth_diameter))
    lon = np.degrees(np.arctan2(y, x))
    return lat,lon

# webscraping(url)
match_link               = fetchData()
data                     = downloadData(match_link)
start_time, finish_time  = headerFinder(data)

for idx in tqdm(range(start_time,finish_time)):
    interpolate_position_x = interpolate(float(data[idx].split()[1]), float(data[idx+1].split()[1]))
    interpolate_position_y = interpolate(float(data[idx].split()[2]), float(data[idx+1].split()[2]))
    interpolate_position_z = interpolate(float(data[idx].split()[3]), float(data[idx+1].split()[3]))
    position_x += [interpolate_position_x]
    position_y += [interpolate_position_y]
    position_z += [interpolate_position_z]
X = np.array(position_x).flatten()
Y = np.array(position_y).flatten()
Z = np.array(position_z).flatten()

for pos in range(len(X)):
    lat,lon = xyz_to_lat_lon(X[pos], Y[pos], Z[pos])
    url = ('https://scihub.copernicus.eu/dhus/search?q=footprint:"Intersects({}, {})"').format(lat, lon)
    print('\n'"ESA Map URL for ISS Current Location: " + url)

# From google
url = ("https://maps.google.com/?q={0},{1}&ll={0},{1}&z=3").format(lat, lon)
print('\n'"Google Map URL for ISS Current Location: " + url)

f=open('static.png','wb')
f.write(requests.get(url).content)
f.close()