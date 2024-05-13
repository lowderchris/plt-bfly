# Script to generate a solar magnetic field butterfly diagram

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import sunpy
import sunpy.io
import sunpy.coordinates
import os
import os.path
import datetime
import pandas

# Grab current sunspot number data from the SILSO dataset
# os.system('curl -O https://www.sidc.be/SILSO/DATA/SN_m_tot_V2.0.txt')
# os.system('curl -O https://www.sidc.be/SILSO/DATA/SN_m_hem_V2.0.txt')

# Define column headers, and read in the data!
dcols = ['yy', 'mm', 't', 'tsn', 'snsd', 'nobs', 'def']
dcolsh = ['yy', 'mm', 't', 'tsn', 'nsn', 'ssn', 'mmsd', 'mmsdn', 'mmsds', 'nobs', 'nobsn', 'nobss', 'prov']

d_tot = pandas.read_csv('SN_m_tot_V2.0.txt', sep='\s+',header=None, names=dcols)
d_hem = pandas.read_csv('SN_m_hem_V2.0.txt', sep='\s+',header=None, names=dcolsh)

# Generate some time arrays
td = []
for i in np.arange(len(d_tot['yy'].values)) : td.append(datetime.datetime(d_tot['yy'].values[i], d_tot['mm'].values[i], 1))
tdh = []
for i in np.arange(len(d_hem['yy'].values)) : tdh.append(datetime.datetime(d_hem['yy'].values[i], d_hem['mm'].values[i], 1))

# Generate an array of corresponding CR numbers
crd = sunpy.coordinates.sun.carrington_rotation_number(td)
crdh = sunpy.coordinates.sun.carrington_rotation_number(tdh)

# Specify data directories and data range
cr0 = 1911
cr1 = 2284
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
        bfly[:,i] = np.interp(lat,np.linspace(-1,1,crdat.shape[0]), 
                              crdat.sum(axis=1)/sum(np.isfinite(crdat[540,:])))
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
tarr_cr = [ datetime.datetime.strptime((tarr[x])[0:16], "%Y.%m.%d_%H:%M") 
           for x in range(len(tarr)) ]
tmin, tmax = mdates.date2num([tarr_cr[0], tarr_cr[-1]])

###

# f, (ax1,ax2) = plt.subplots(2,1, figsize=[6,6], gridspec_kw=kw, sharex=True)

# df.plot(kind='scatter', x='x',  y='y', c='score', s=100, cmap="PuRd",
#           ax=ax, colorbar=True)
# df.groupby("x").mean().plot(kind = 'bar', y='score',ax=ax2, legend=False)

# ax2.legend(bbox_to_anchor=(1.03,0),loc=3)

# pos = ax.get_position()
# pos2 = ax2.get_position()
# ax2.set_position([pos.x0,pos2.y0,pos.width,pos2.height])

##

f, (ax1, ax2) = plt.subplots(2, figsize=[6,5], gridspec_kw={'height_ratios':[2,1]}, sharex=True)
im = ax1.imshow(bfly, vmin=-10, vmax=10, extent=[tmin,tmax,-1,1], 
                aspect='auto', cmap='Greys_r')
ax1.set_yticks([-1,-0.5,0,0.5,1])
ax1a = ax1.twinx()
ax1b = ax1.twiny()
latticks = np.array([-90,-60,-30,30,60,90])
ax1a.set_yticks(np.sin(latticks*np.pi/180))
ax1a.set_yticklabels(latticks)
ax1b.set_xlim(cr0, cr1)
ax1.xaxis_date()
ax1.set_ylabel('Sine latitude')
ax1a.set_ylabel('Latitude (degrees)')
ax1b.set_xlabel('Carrington rotation')
cb = plt.colorbar(im, ax=ax1a, label='Mean magnetic flux density (G)', 
                  use_gridspec=True, fraction=0.1, pad=0.15, extend='both')

ax2.plot(mdates.date2num(td), d_tot['tsn'], 'k', label='Total')
ax2.plot(mdates.date2num(tdh), d_hem['nsn'], 'C0', label='Northern')
ax2.plot(mdates.date2num(tdh), d_hem['ssn'], 'C3', label='Southern')
ax2.set_xlim(tmin,tmax)
ax2b = ax2.twiny()
ax2b.set_xlim(cr0, cr1)
ax2b.set_xticklabels([])
ax2.set_ylim(0,275)
ax2.legend(bbox_to_anchor=(1.03,0),loc=3)
ax2.xaxis_date()
ax2.set_xlabel('Year')
ax2.set_ylabel('Sunspot number')

# vline_crs = [2130, 2160, 2193, 2231]
vline_crs = []
ax1b.vlines(vline_crs, -1, 1, linestyles='dashed', color='k', alpha=0.5)
ax2b.vlines(vline_crs, 0, 275, linestyles='dashed', color='k', alpha=0.5)

pos = ax1.get_position()
pos2 = ax2.get_position()
ax2.set_position([pos.x0,pos2.y0,pos.width,pos2.height])

#f.tight_layout()
plt.savefig('bfly_ssn.pdf', dpi=300)
plt.savefig('bfly_ssn.png', dpi=300)

##

f, (ax1, ax2) = plt.subplots(2, figsize=[9,5], gridspec_kw={'height_ratios':[2,1]}, sharex=True)
im = ax1.imshow(bfly, vmin=-10, vmax=10, extent=[tmin,tmax,-1,1], 
                aspect='auto', cmap='Greys_r')
ax1.set_yticks([-1,-0.5,0,0.5,1])
ax1a = ax1.twinx()
ax1b = ax1.twiny()
latticks = np.array([-90,-60,-30,30,60,90])
ax1a.set_yticks(np.sin(latticks*np.pi/180))
ax1a.set_yticklabels(latticks)
ax1b.set_xlim(cr0, cr1)
ax1.xaxis_date()
ax1.set_ylabel('Sine latitude')
ax1a.set_ylabel('Latitude (degrees)')
ax1b.set_xlabel('Carrington rotation')
cb = plt.colorbar(im, ax=ax1a, label='Mean magnetic flux density (G)', 
                  use_gridspec=True, fraction=0.05, pad=0.1, extend='both')

ax2.plot(mdates.date2num(td), d_tot['tsn'], 'k', label='Total')
ax2.plot(mdates.date2num(tdh), d_hem['nsn'], 'C0', label='Northern')
ax2.plot(mdates.date2num(tdh), d_hem['ssn'], 'C3', label='Southern')
ax2.set_xlim(tmin,tmax)
ax2b = ax2.twiny()
ax2b.set_xlim(cr0, cr1)
ax2b.set_xticklabels([])
ax2.set_ylim(0,275)
ax2.legend(bbox_to_anchor=(1.03,0),loc=3)
ax2.xaxis_date()
ax2.set_xlabel('Year')
ax2.set_ylabel('Sunspot number')

# vline_crs = [2130, 2160, 2193, 2231]
vline_crs = []
ax1b.vlines(vline_crs, -1, 1, linestyles='dashed', color='k', alpha=0.5)
ax2b.vlines(vline_crs, 0, 275, linestyles='dashed', color='k', alpha=0.5)

pos = ax1.get_position()
pos2 = ax2.get_position()
ax2.set_position([pos.x0,pos2.y0,pos.width,pos2.height])

#f.tight_layout()
plt.savefig('bfly_ssn2.pdf', dpi=300)
plt.savefig('bfly_ssn2.png', dpi=300)