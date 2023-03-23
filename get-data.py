# Script to grab and update MDI and HMI data
# See sunpy notes at https://sunpy.readthedocs.io/en/latest/guide/acquiring_data/jsoc.html

# Import libraries
import sunpy
import sunpy.io
import sunpy.coordinates
from sunpy.net import Fido, attrs as a
import drms
import os
import glob

# Specify any directories
hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil/')
mdidat = os.path.expanduser('~/data/mdi.Synoptic_Mr.polfil/')

# HMI data

# Sort out the last downloaded rotation
crfiles = glob.glob(hmidat+'*.fits')
crfiles.sort()
crlist = [int(i[-19:-15]) for i in crfiles]

# Specify requested rotations
if len(crlist) == 0:
    cr0 = 2096
else:
    cr0 = max(crlist) + 1
cr1 = int(sunpy.coordinates.sun.carrington_rotation_number(t='now'))

if (cr0 - 1) == cr1:
    print('what?')

# Start the client
c = drms.Client()

# Generate a search
crots = a.jsoc.PrimeKey('CAR_ROT', str(cr0) + '-' + str(cr1))
res = Fido.search(a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'), crots, 
                  a.jsoc.Notify(os.environ["JSOC_EMAIL"]))

# Once the query is made and trimmed down...
download = Fido.fetch(res, path=hmidat+'{file}.fits')

# MDI data

# Grab MDI
os.system('mkdir ' + mdidat)
os.chdir(mdidat)
os.system('curl -O "http://soi.stanford.edu/magnetic/synoptic/carrot/M_Corr/synop_Mr_0.polfil.[1911-2104].fits"')