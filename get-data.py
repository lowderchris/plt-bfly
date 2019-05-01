# Script to grab and update MDI and HMI data
# See sunpy notes at https://sunpy.readthedocs.io/en/latest/guide/acquiring_data/jsoc.html

# Import libraries
import sunpy
import datetime
from sunpy.net import Fido, attrs as a
import drms

# Specify any directories
hmidat = '/Users/clowder/data/hmi.Synoptic_Mr.polfil/'

# Start the client
c = drms.Client()

# Generate a search
# CL - sanitize the use of my email address and directories here...
today = datetime.datetime.now().replace(microsecond=0,second=0,minute=0,hour=0)
res = Fido.search(a.Time('2009-01-01', today), a.jsoc.Series('hmi.Synoptic_Mr_polfil_720s'), a.jsoc.Notify(''))


# Once the query is made and trimmed down...
download = Fido.fetch(res, path=hmidat+'{file}.fits')
