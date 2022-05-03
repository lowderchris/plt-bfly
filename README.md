# plt-bfly
Generate a solar magnetic field butterfly diagram

![alt text](https://github.com/lowderchris/plt-bfly/blob/master/bfly.png?raw=true)

## Grabbing data
HMI data can be acquired using JSOC, or other similar tools. The included script will automatically fetch HMI polar corrected synoptic charts, with a few specifications of directory structure, JSOC email, etc.

To access the set of polar corrected MDI synoptic charts:

    curl -O http://soi.stanford.edu/magnetic/synoptic/carrot/M_Corr/synop_Mr_0.polfil.[1911-2104].fits
