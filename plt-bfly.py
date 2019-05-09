# Script to generate a solar magnetic field butterfly diagram

# Import libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sunpy
import sunpy.io
import os.path

# Specify data directories and data range
mdidat = os.path.expanduser('~/data/mdi.Synoptic_Mr.polfil/')
mdicrs = np.arange(1911, 2096)
hmidat = os.path.expanduser('~/data/hmi.Synoptic_Mr.polfil/')
hmicrs = np.arange(2097, 2215)

# Specify any parameters needed for the output map
latres = 1440
lat = np.linspace(-1,1,latres)

# Generate any output arrays
bfly = np.zeros([latres, len(mdicrs)+len(hmicrs)], dtype=np.double)
i = 0

# Loop through MDI files
for cr in mdicrs:
    # Read and process a given synoptic chart into a profile
    filenm = mdidat+'synop_Mr_0.polfil.'+str(cr).zfill(4)+'.fits'
    if os.path.isfile(filenm):
        d = sunpy.io.read_file(filenm)
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
        crdat = d[1].data
        bfly[:,i] = crdat.sum(axis=1)/crdat.shape[1]
    else:
        bfly[:,i] = np.NaN
    i += 1

# Plot the resulting diagram
f, (ax1) = plt.subplots(1, figsize=[8,3])
im = ax1.imshow(bfly, vmin=-10, vmax=10, extent=[mdicrs[0],hmicrs[-1],-1,1], aspect='auto')
ax1.set_yticks([-1,-0.5,0,0.5,1])
ax2 = ax1.twinx()
latticks = np.array([-90,-60,-30,30,60,90])
ax2.set_yticks(np.sin(latticks*np.pi/180))
ax2.set_yticklabels(latticks)
ax1.set_xlabel('Carrington rotation')
ax1.set_ylabel('Sine latitude')
ax2.set_ylabel('Latitude (degrees)')
cb = plt.colorbar(im, ax=ax2, label='Mean magnetic field (G)', use_gridspec=True, fraction=0.05, pad=0.1, extend='both')
f.tight_layout()
plt.savefig('bfly.pdf', dpi=300)
plt.savefig('bfly.png', dpi=300)
