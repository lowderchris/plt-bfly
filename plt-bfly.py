# Script to generate a solar magnetic field butterfly diagram

# Import libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sunpy
import sunpy.io
import os.path
import datetime

# Specify data directories and data range
cr0 = 1911
cr1 = 2231
mdidat = os.path.expanduser('~/data/mdi.Synoptic_Mr.polfil/')
mdicrs = np.arange(cr0, 2096)
hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil/')
hmicrs = np.arange(2097, cr1)

# Specify any parameters needed for the output map
latres = 1440
lat = np.linspace(-1,1,latres)

# Generate a time array
tarr = []

# Generate any output arrays
bfly = np.zeros([latres, len(mdicrs)+len(hmicrs)], dtype=np.double)
i = 0

# Loop through MDI files
for cr in mdicrs:
    # Read and process a given synoptic chart into a profile
    filenm = mdidat+'synop_Mr_0.polfil.'+str(cr).zfill(4)+'.fits'
    if os.path.isfile(filenm):
        d = sunpy.io.read_file(filenm)
        tarr.append(d[0].header['T_ROT'])
        crdat = d[0].data
        crdat[np.where(abs(crdat) > 3000)] = 0
        bfly[:,i] = np.interp(lat,np.linspace(-1,1,crdat.shape[0]),crdat.sum(axis=1)/sum(np.isfinite(crdat[540,:])))
    else:
        bfly[:,i] = np.NaN
    i += 1

# Loop through HMI files
for cr in hmicrs:
    # Read and process a given synoptic chart into a profile
    filenm = hmidat+'hmi.synoptic_mr_polfil_720s.'+str(cr).zfill(4)+'.Mr_polfil.fits'
    if os.path.isfile(filenm):
        d = sunpy.io.read_file(filenm)
        tarr.append(d[1].header['T_ROT'])
        crdat = d[1].data
        bfly[:,i] = crdat.sum(axis=1)/crdat.shape[1]
    else:
        bfly[:,i] = np.NaN
    i += 1

# Convert the timing array
tarr_cr = [ datetime.datetime.strptime((tarr[x])[0:16], "%Y.%m.%d_%H:%M") for x in range(len(tarr)) ]
tmin, tmax = mdates.date2num([tarr_cr[0], tarr_cr[-1]])

# Plot the resulting diagram
f, (ax1) = plt.subplots(1, figsize=[8,4])
im = ax1.imshow(bfly, vmin=-10, vmax=10, extent=[tmin,tmax,-1,1], aspect='auto', cmap='Greys_r')
ax1.set_yticks([-1,-0.5,0,0.5,1])
ax1a = ax1.twinx()
ax1b = ax1.twiny()
latticks = np.array([-90,-60,-30,30,60,90])
ax1a.set_yticks(np.sin(latticks*np.pi/180))
ax1a.set_yticklabels(latticks)
ax1b.set_xlim(cr0, cr1)
ax1.xaxis_date()
ax1.set_xlabel('Year')
ax1.set_ylabel('Sine latitude')
ax1a.set_ylabel('Latitude (degrees)')
ax1b.set_xlabel('Carrington rotation')
cb = plt.colorbar(im, ax=ax1a, label='Mean magnetic field (G)', use_gridspec=True, fraction=0.05, pad=0.1, extend='both')
#f.tight_layout()
plt.savefig('bfly.pdf', dpi=300)
plt.savefig('bfly.png', dpi=300)
