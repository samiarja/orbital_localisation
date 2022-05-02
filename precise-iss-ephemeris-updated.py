from typing import Final
import httplib2
import requests
import os
from time import time, sleep
import numpy as np
from tqdm import tqdm
from datetime import datetime, timezone
from bs4 import BeautifulSoup, SoupStrainer
import urllib.request
from PIL import Image
import math


'''
https://scihub.copernicus.eu/twiki/do/view/SciHubUserGuide/FullTextSearch?redirectedfrom=SciHubUserGuide.3FullTextSearch
https://scihub.copernicus.eu/dhus/#/home
'''

position_x      = []
position_y      = []
position_z      = []
earth_diameter  = 6371
ZOOM            = 6
START  = '2022-01-24T20:07:00'
FINISH = '2022-01-24T20:15:00'
url = 'https://data.nasa.gov/browse?q=ISS+COORDS&sortBy=newest&utf8=%E2%9C%93&page='

class GoogleMapsLayers:
  ROADMAP = "v"
  TERRAIN = "p"
  ALTERED_ROADMAP = "r"
  SATELLITE = "s"
  TERRAIN_ONLY = "t"
  HYBRID = "y"


class GoogleMapDownloader:
    """
        A class which generates high resolution google maps images given
        a longitude, latitude and zoom level
    """

    def __init__(self, lat, lng, zoom=12, layer=GoogleMapsLayers.ROADMAP):
        """
            GoogleMapDownloader Constructor
            Args:
                lat:    The latitude of the location required
                lng:    The longitude of the location required
                zoom:   The zoom level of the location required, ranges from 0 - 23
                        defaults to 12
        """
        self._lat = lat
        self._lng = lng
        self._zoom = zoom
        self._layer = layer

    def getXY(self):
        """
            Generates an X,Y tile coordinate based on the latitude, longitude
            and zoom level
            Returns:    An X,Y tile coordinate
        """

        tile_size = 256

        # Use a left shift to get the power of 2
        # i.e. a zoom level of 2 will have 2^2 = 4 tiles
        numTiles = 1 << self._zoom

        # Find the x_point given the longitude
        point_x = (tile_size / 2 + self._lng * tile_size / 360.0) * numTiles // tile_size

        # Convert the latitude to radians and take the sine
        sin_y = math.sin(self._lat * (math.pi / 180.0))

        # Calulate the y coorindate
        point_y = ((tile_size / 2) + 0.5 * math.log((1 + sin_y) / (1 - sin_y)) * -(
        tile_size / (2 * math.pi))) * numTiles // tile_size

        return int(point_x), int(point_y)

    def generateImage(self, **kwargs):
        """
            Generates an image by stitching a number of google map tiles together.

            Args:
                start_x:        The top-left x-tile coordinate
                start_y:        The top-left y-tile coordinate
                tile_width:     The number of tiles wide the image should be -
                                defaults to 5
                tile_height:    The number of tiles high the image should be -
                                defaults to 5
            Returns:
                A high-resolution Goole Map image.
        """

        start_x = kwargs.get('start_x', None)
        start_y = kwargs.get('start_y', None)
        tile_width = kwargs.get('tile_width', 5)
        tile_height = kwargs.get('tile_height', 5)

        # Check that we have x and y tile coordinates
        if start_x == None or start_y == None:
            start_x, start_y = self.getXY()

        # Determine the size of the image
        width, height = 256 * tile_width, 256 * tile_height

        # Create a new image of the size require
        map_img = Image.new('RGB', (width, height))

        for x in range(0, tile_width):
            for y in range(0, tile_height):
                url = f'https://mt0.google.com/vt?lyrs={self._layer}&x=' + str(start_x + x) + '&y=' + str(start_y + y) + '&z=' + str(
                    self._zoom)

                current_tile = str(x) + '-' + str(y)
                urllib.request.urlretrieve(url, current_tile)

                im = Image.open(current_tile)
                map_img.paste(im, (x * 256, y * 256))

                os.remove(current_tile)

        return map_img

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
    position_interpolation = np.linspace(float(x1), float(x2), num=20)
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

for pos in tqdm(range(len(X))):
    lat,lon = xyz_to_lat_lon(X[pos], Y[pos], Z[pos])

    # gmd = GoogleMapDownloader(lat,lon, ZOOM, GoogleMapsLayers.SATELLITE)
    gmd = GoogleMapDownloader(lat,lon, ZOOM, GoogleMapsLayers.HYBRID)
    print("The tile coordinates are {}".format(gmd.getXY()))

    try:
        # Get the high resolution image
        img = gmd.generateImage()
    except IOError:
        print("Could not generate the image - try adjusting the zoom level and checking your coordinates")
    else:
        # Save the image to disk
        img.save("MAPS/map_HYBRID" + str(pos) + ".png")
        print("The map has successfully been created")

'''
Convert 
ffmpeg -r 10 -f image2 -s 1920x1080 -i map_HYBRID%01d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p mapGenerationVideo.mp4
'''
