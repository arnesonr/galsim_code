import galsim
import csv
import numpy
import math

chifilename = '/home/rarneson/new_GalSim/gal_R22_S11.dat'
a_im = []
b_im = []
ra_im = []
dec_im =[]
flux_im = []
theta_im = []
sersic_im = []

with open(chifilename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader:
        a_im.append(float(row[6]))
        b_im.append(float(row[7]))
        ra_im.append(float(row[2]))
        dec_im.append(float(row[3]))
        flux_im.append(float(row[5]))
        theta_im.append(float(row[8]))
        sersic_im.append(float(row[9]))
#From Chihway's config file
photon=0
Pixelsize=0.2
full_image = galsim.ImageF(128, 128)
full_image.setCenter(2657,3044)
full_image.setScale(0.2)
im_center = full_image.bounds.center()
im_center = galsim.PositionD(im_center.x, im_center.y)
Atmseeing = 0.6
Airmass = 1.0
atm_psf = galsim.Kolmogorov(fwhm = Atmseeing*Airmass**0.6)
Optpsfbeta = 3
Optpsfsize = 0.35
Optpsftrunc = 3.*Optpsfsize
opt_psf=galsim.Moffat(beta=Optpsfbeta, fwhm=Optpsfsize, trunc=Optpsftrunc)
pix = galsim.Pixel(0.2)
#decide if using photon shooting or not
psf=galsim.Convolve([atm_psf, opt_psf, pix])    


if (photon == 0):
    psf=galsim.Convolve([atm_psf, opt_psf, pix])
if (photon == 1):
    psf=galsim.Convolve([atm_psf, opt_psf])

#Bright galaxy
#corresponds to line 24645 in storedproc.dat
x=512.84819306
y=652.163783834
z=0.215078294
galFlux=int(19861.5428849*25)
galA=3.0203414
galB=3.01904535
galPhi=118.723953+numpy.pi/2
# this will change in the next phosim version!!
galN=4.0
re=(galA*galB)**0.5
gale=(galA**2-galB**2)/(galA**2+galB**2)

#Another galaxy
#corresponds to 14799 in storedproc.dat

x=2656.87679306
y=3044.47142383
z=0.0981191993
galFlux=int(154.754469487*25)
galA=0.222630903
galB=0.169356003
galPhi=309.135712+numpy.pi/2
# this will change in the next phosim version!!
galN=4.0
re=(galA*galB)**0.5
gale=(galA**2-galB**2)/(galA**2+galB**2)

x=1816.68989306
y=2161.70692783
z=1.24828625
galFlux=int(1.95640690433*25)
galA=0.496319503
galB=0.471401006
galPhi=355.426819+numpy.pi/2
# this will change in the next phosim version!!
galN=4.0
re=(galA*galB)**0.5
gale=(galA**2-galB**2)/(galA**2+galB**2)

# Make the galaxy profile with these values:
gal = galsim.Sersic(half_light_radius=re, flux=galFlux, n=galN, gsparams = galsim.GSParams(maximum_fft_size = 10000))
gal.applyShear(e=gale, beta=galPhi*galsim.radians)

# Get the integer values of these which will be the center of the 
# postage stamp image.
ix = int(math.floor(x+0.5))
iy = int(math.floor(y+0.5))
# The remainder will be accounted for in a shift
dx = x - ix
dy = y - iy
# Build the final object
final = galsim.Convolve([psf, gal])
# Account for the non-integral portion of the position
final.applyShift(dx*Pixelsize,dy*Pixelsize)

# Draw the stamp image
if (photon==0):
stamp = final.draw(wmult=5, dx=Pixelsize)
if (photon==1):
stamp=final.drawShoot(wmult=5, n_photons=galFlux, dx=Pixelsize)

# Recenter the stamp at the desired position:
stamp.setCenter(ix,iy)
# Find overlapping bounds
bounds = stamp.bounds & full_image.bounds
full_image[bounds] += stamp[bounds]

x=512.84819306
y=652.163783834
z=0.215078294
galFlux=int(2153.71072635*25)
galA=5.11530113
galB=1.96511185
galPhi=118.723953+numpy.pi/2
# this will change in the next phosim version!!
galN=1.0
re=(galA*galB)**0.5
gale=(galA**2-galB**2)/(galA**2+galB**2)

#corresponds to line 14799 in storedproc.dat (use '-x 29148. -y 4382. --only-line 14798'in galsimcat.py)
x=1816.68989306
y=2161.70692783
z=1.24828625
galFlux=int(1324.55115068*25)
galA=0.587640882
galB=0.143685102
galPhi=355.426819+numpy.pi/2
# this will change in the next phosim version!!
galN=1.0
re=(galA*galB)**0.5
gale=(galA**2-galB**2)/(galA**2+galB**2)

#corresponds to line 13 in storedproc.dat (use '-x 28308 -y 3499 --only-line 13403'in galsimcat.py)
x=2656.87679306
y=3044.47142383
z=0.0981191993
galFlux=int(10277.4074546*25)
galA=5.24951363
galB=1.24920821
galPhi=309.135712+numpy.pi/2
# this will change in the next phosim version!!
galN=1.0
re=(galA*galB)**0.5
gale=(galA**2-galB**2)/(galA**2+galB**2)

# Make the galaxy profile with these values:
gal = galsim.Sersic(half_light_radius=re, flux=galFlux, n=galN, gsparams = galsim.GSParams(maximum_fft_size = 10000))
gal.applyShear(e=gale, beta=galPhi*galsim.radians)

# Get the integer values of these which will be the center of the 
# postage stamp image.
ix = int(math.floor(x+0.5))
iy = int(math.floor(y+0.5))
# The remainder will be accounted for in a shift
dx = x - ix
dy = y - iy
# Build the final object
final = galsim.Convolve([psf, gal])
# Account for the non-integral portion of the position
final.applyShift(dx*Pixelsize,dy*Pixelsize)

# Draw the stamp image
if (photon==0):
stamp = final.draw(wmult=5, dx=Pixelsize)
if (photon==1):
stamp=final.drawShoot(wmult=5, n_photons=galFlux, dx=Pixelsize)

# Recenter the stamp at the desired position:
stamp.setCenter(ix,iy)
# Find overlapping bounds
bounds = stamp.bounds & full_image.bounds
full_image[bounds] += stamp[bounds]

rng = galsim.UniformDeviate(145604320) # Initialize the random number generator 
full_image.addNoise(galsim.PoissonNoise(rng, sky_level=160.901631573))
full_image.write('chihway_gal.fits', clobber = True)


imsim = galsim.fits.read('R22_S11.fits')
gsim = galsim.fits.read('storedproc_noise.fits')
diff_im = imsim-gsim
diff_im.write('diff_im.fits', clobber = True)

###
#GalSim
####
filename='/home/rarneson/new_GalSim/storedproc.dat'
ra, dec, pa_b, pa_d, a_b, a_d, b_b, b_d, hlr_b, hlr_d,flux_b,flux_d = [],[],[],[],[],[],[],[],[],[],[],[]
flux_row = 23
mean_lambda = 754.5 
dlambda = 127.
with open(filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    reader.next()
    for row in reader:
        ra.append(float(row[1]))
        dec.append(float(row[2]))
        pa_b.append(float(row[7]))
        pa_d.append(float(row[8]))
        a_b.append(float(row[15]))
        a_d.append(float(row[16]))
        b_b.append(float(row[17]))
        b_d.append(float(row[18]))
        hlr_b.append(float(row[5]))
        hlr_d.append(float(row[6]))
        flux_b.append(float(row[9]))
        flux_d.append(float(row[10]))

#compare Bluge component
for i in range(n.size(a_b)):
    for j in range(n.size(a_im)):
        if n.abs(a_b[i]-a_im[j]) < 0.00000001 and n.abs(b_b[i]-b_im[j]) < 0.00000001:
            print i,j
#compare Disk component
for i in range(n.size(a_d)):
    for j in range(n.size(a_im)):
        if n.abs(a_d[i]-a_im[j]) < 0.00000001 and n.abs(b_d[i]-b_im[j]) < 0.00000001:
            print i,j
me =[]
me_ra=[]
chihway =[]
chihway_ra=[]
samefilename='/home/rarneson/new_GalSim/samelistd.txt'
with open(samefilename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    for row in reader:
        me.append(int(row[0]))
        chihway.append(int(row[1]))

n.size(me)
n.size(chihway)
how_many = n.size(me)
for i in range(0,how_many):
    me_ra.append(ra[me[i]])
    chihway_ra.append(ra_im[chihway[i]])

new_me_ra=[]
for i in range(n.size(me_ra)):
    new_me_ra.append((me_ra[i] - 1.58285320483)*18000.+2036.)

