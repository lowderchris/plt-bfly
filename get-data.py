# Script to grab and update MDI and HMI data
# See sunpy notes at https://sunpy.readthedocs.io/en/latest/guide/acquiring_data/jsoc.html

# Import libraries
import sunpy
import sunpy.io
import datetime
from sunpy.net import Fido, attrs as a
import drms
import os

# Specify any directories
hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil/')

# Start the client
c = drms.Client()

# Generate a search
# CL - sanitize the use of my email address and directories here...
today = datetime.datetime.now().replace(microsecond=0,second=0,minute=0,hour=0)
res = Fido.search(a.Time('2019-06-01', today), a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'))


# Once the query is made and trimmed down...
download = Fido.fetch(res, path=hmidat+'{file}.fits')
