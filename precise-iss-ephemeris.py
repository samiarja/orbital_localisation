from IPython.display import HTML
import numpy as np
from matplotlib import animation, gridspec
import matplotlib.pyplot as plt
from datetime import datetime, timezone
import requests
import urllib.request
import os
from requests import head
import time
import scipy.io as sio
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D

'''
ISS COORDS Archive
https://data.nasa.gov/browse?q=ISS%20COORDS&sortBy=relevance
'''

'FROM 2022-01-24T20:11:00.000 TO 2022-01-24T20:39:00.000'


position_x      = []
position_y      = []
position_z      = []
trajectory_data = "ISS_TOPO/ISS.OEM_J2K_EPH.txt"
image_file = 'img/blue_marble.jpg'
url = "https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.txt"

if not os.path.isdir("ISS_TOPO"):
    os.mkdir("ISS_TOPO")

r = requests.get(url, allow_redirects=True)
open('ISS_TOPO/ISS.OEM_J2K_EPH.txt', 'wb').write(r.content)

with open(trajectory_data) as file:
    lines = [line.rstrip() for line in file]

def connection(url='https://spotthestation.nasa.gov/trajectory_data.cfm', timeout=5):
    try:
        req = requests.get(url, timeout=timeout)
        req.raise_for_status()
        return True
    except requests.HTTPError as e:
        print("Checking internet connection failed, status code {0}.".format(
        e.response.status_code))
    except requests.ConnectionError:
        print("No internet connection available.")
    return False

def headerFinder(infile):
    for num,line in enumerate(infile):
        if line.startswith("COMMENT End sequence of events"):
            return num

def mpl_sphere(image_file):
    R = 5000
    img = plt.imread(image_file)

    theta = np.linspace(0, np.pi, img.shape[0])
    phi = np.linspace(0, 2*np.pi, img.shape[1])

    count = 180
    theta_inds = np.linspace(0, img.shape[0] - 1, count).round().astype(int)
    phi_inds = np.linspace(0, img.shape[1] - 1, count).round().astype(int)
    theta = theta[theta_inds]
    phi = phi[phi_inds]
    img = img[np.ix_(theta_inds, phi_inds)]
    theta,phi = np.meshgrid(theta, phi)
    x = R * np.sin(theta) * np.cos(phi)
    y = R * np.sin(theta) * np.sin(phi)
    z = R * np.cos(theta)
    return x, y, z, img

def interpolate(x1,x2):
    position_interpolation = np.linspace(float(x1), float(x2), num=4)
    return position_interpolation

def visualise_glob(img,glob_x,glob_y,glob_z):
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(glob_x.T, glob_y.T, glob_z.T, facecolors=img/255, cstride=1, rstride=1)
    ax.scatter(position_x[:], position_y[:], position_z[:], alpha=0.2)
    ax.set_title("ISS Trajectory Path", fontsize=16)
    ax.set_xlabel("Position x", fontsize=16)
    ax.set_ylabel("Position y", fontsize=16)
    ax.set_zlabel("Position z", fontsize=16)
    ax.axis('auto')
    plt.tight_layout()
    plt.show()

def main(plot=False):
    header = headerFinder(lines)
    for idx in tqdm(range(header+1,len(lines)-1)):
        first_line = lines[idx]
        lines_split = first_line.split()
        timeStamp = lines_split[0]
        timeStamp_formatted = timeStamp[:-4]
        date_conversion = datetime.strptime(timeStamp_formatted, "%Y-%m-%dT%H:%M:%S")
        interpolate_position_x = interpolate(float(lines[idx].split()[1]), float(lines[idx+1].split()[1]))
        interpolate_position_y = interpolate(float(lines[idx].split()[2]), float(lines[idx+1].split()[2]))
        interpolate_position_z = interpolate(float(lines[idx].split()[3]), float(lines[idx+1].split()[3]))
        position_x.append(interpolate_position_x)
        position_y.append(interpolate_position_y)
        position_z.append(interpolate_position_z)
    glob_x, glob_y, glob_z, img = mpl_sphere(image_file)
    if plot:
        visualise_glob(img,glob_x,glob_y,glob_z)

main(plot=True)
time.sleep(60)

if connection() == True:
    print("ISS TOPO Tracker (OEM)")
    try:
        while True:
            main()
    except(KeyboardInterrupt):
        print("Programme Interrupted")
else:
    quit()