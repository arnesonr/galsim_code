import os
import numpy as n
import pyfits as py
from numpy import matrix
from numpy import linalg

def neff(line, inname, outname, flux_norm = 711.0, skylevel=1020.67, pixel_scale=0.2, width='4096', height='4096', x='0.5', y='0.0', N_exp=1.):
	"""
	ARGUMENTS:
      line: which line to send to galsimcat *NOTE: using 0 for line will draw all galaxies in region 
      inname: [string] the input dat filename to send to galsimcat
      outname: [string] the output filename
      flux_norm: total flux in ADU for one exposure of a typical galaxy of AB mag 24
      skylevel: the skylevel to send to galsimcat and used in calculating neff, default is i-band
                NOTE: average sky level is (ADU/pixel in ONE exposure)
      pixel_scale: the size of each pixel in arcsec
      height/width: pixel size of image (0.2 arcsec/pixel)
      x/y: center of image in degrees
      N_exp: number of exposures (i & r band will have 460 (15 sec. each) exposures in 10 yr. span)
	"""
	#"""
	os.system('./galsimcat.py --stamps --partials -i '+str(inname)+' -o' +str(outname)+' --only-line '
		+str(line)+' --flux-norm '+str(flux_norm)+' --sky-level '+str(skylevel)+' --pixel-scale '+str(pixel_scale)+' --width '+str(width)+' --height '+str(height)+
		' --margin 5.0 -x '+str(x)+' -y '+str(y)+' --no-trim')
	#"""
	#first entry in 7x7 (or 6x6 if no bulge) matrix is flux which we
	#can get from nominal
	n_eff=0.0
	hdulist = py.open(outname+'_stamps.fits',memmap = False)
	N_gal = n.shape(hdulist)[0]
	#make a txt file to write neff to
	os.system('touch '+str(outname)+'.txt')
	f = open(str(outname)+'.txt','w')
	for k in range(1,N_gal):
		scidata = hdulist[k].data
		nominal = scidata[0].flatten()
		xc = scidata[3].flatten()
		yc = scidata[4].flatten()
		hlr_d = scidata[5].flatten()
		hlr_b = scidata[6].flatten()
		g1 = scidata[7].flatten()
		g2 = scidata[8].flatten()
		cube = list()
		cube.append(nominal)
		cube.append(xc)
		cube.append(yc)
		if n.sum(hlr_d) != 0.0 and n.sum(hlr_b) != 0.0:
			cube.append(hlr_d)
			cube.append(hlr_b)
			mat_size = 7
			MAT = matrix(n.zeros((7,7)))
			#MAT = matrix([[0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.,0.]])

		if n.sum(hlr_d) == 0.0 and n.sum(hlr_b) != 0.0:
			cube.append(hlr_b)
			mat_size = 6
			MAT = matrix(n.zeros((6,6)))
			#MAT = matrix([[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.,0.]])

		if n.sum(hlr_d) != 0.0 and n.sum(hlr_b) == 0.0:
			cube.append(hlr_d)
			mat_size = 6
			MAT = matrix(n.zeros((6,6)))
			
		cube.append(g1)
		cube.append(g2)
		#need to calculate the parameter crosscorelation
		#the variance in each pixel is equal to the sky_level + pixel flux
		sigma_SN = 0.25
		holder = n.copy(nominal)*0.0
		for i in range(mat_size):
			for j in range(mat_size):
				MAT[i,j] = n.sum((cube[i]*cube[j])/(skylevel+cube[0]))
		
		if linalg.det(MAT) != 0.0:
			#invert the matrix
			IMAT = MAT.I
		if linalg.det(MAT) == 0.0:
			print 'Singular matrix encountered skipping line:', k
			f.write('0\n')
			continue
		#find the determinant of the sub matrix {e1,e2}
		SUBM = matrix([[IMAT[(mat_size-2),(mat_size-2)],IMAT[(mat_size-2),(mat_size-1)]],[IMAT[(mat_size-1),(mat_size-2)],IMAT[(mat_size-1),(mat_size-1)]]])
		variance_m = (SUBM[0,0] + SUBM[1,1])/2.0
		#At this stage I wont bother with writing out the galaxy id, 
		#as it can be found in the matching ..._catalog.dat file
		#f.write(str(SUBM[0,0])+','+str(SUBM[1,1])+'\n')
		f.write(str(variance_m)+'\n') # python will convert \n to os.linesep
		
		#print 'neff(',k,'):', sigma_SN**2/(sigma_SN**2 + variance_m)
		
		n_eff += sigma_SN**2/(sigma_SN**2 + (variance_m/N_exp))
	#on the last line print the total neff
	f.write(str(n_eff)+'\n')
	f.close()
	hdulist.close()
	return n_eff

