# Script to generate a solar magnetic field butterfly diagram

# Import libraries
import sunpy
import sunpy.io

# Specify data directories and data range
hmidat = '/Users/clowder/data/hmi.Synoptic_MR_720s/'
crs = np.arange(2097, 2213)

# Generate any output arrays
bfly = np.zeros([1440, len(crs)], dtype=np.double)
i = 0

# Loop
for cr in crs:

    # Read and process a given synoptic chart into a profile
    d = sunpy.io.read_file(hmidat+'hmi.Synoptic_Mr_720s.'+str(cr).zfill(4)+'.synopMr.fits')
    crdat = d[0].data
    bfly[:,i] = crdat.sum(axis=1)/3600
    i += 1
