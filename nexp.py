import numpy as n
import pyfits as py
import galsim
from numpy import matrix
import galsimcat
import math

def exposure_plot(inname, sigma_SN = 0.25):
    f = open(str(inname)+'.txt','r')
    variance_m = []
    for line in f:
        if line != 0:
            variance_m.append(float(line))
    #take off the last entry in inname as it is the neff
    variance_m.pop()
    N_exp = n.arange(1,921)*1.0
    n_eff = n.copy(N_exp)*0.0
    for j in range(n.size(N_exp)):
        this_neff = 0.0
        for i in range(n.size(variance_m)):
            this_neff += sigma_SN**2/(sigma_SN**2 + (variance_m[i]/N_exp[j]))
        n_eff[j] = this_neff
    return n_eff

def optimal_weight(inname, skylevel = 1020.67, sigma_SN = 0.25, N_exp = 1, size=1024):
    #find and store all the bbox
    f = open(str(inname)+'_catalog.dat','r')
    bounds= list()
    lines=[]
    for line in f:
        lines.append(line.split())
    f.close
    how_many = n.size(lines)/n.size(lines[0])
    #extract all the bounds
    for i in range(how_many):
        this_bounds = galsim.BoundsI(int(lines[i][8]),int(lines[i][9]),int(lines[i][10]),int(lines[i][11]))
        bounds.append(this_bounds)
    #loop through galaxies and see which overlap
    matches = list()
    for i in range((how_many)):
        for j in range(how_many):
            if ((bounds[i] & bounds[j]).area()) > 0.0:
                matches.append((i,j))
                #these are the only sections of the matrix we need to evaluate
                #all other parts are zero
                #even if there are no overlaps matches will account for diagonal elements
    
    MAT = matrix(n.zeros(((7*how_many),(7*how_many))))

    #we need to go through all the stamps and calcualte matrix elements
    #we need the stamp positions relative to each other
    #loop over all the matches
    hdulist = py.open(str(inname)+'_stamps.fits',memmap=False)
    #make a list of the elements to delete
    deletelist = list()
    tokeepornot = list()
    n_effi = 0.0
    n_effc = 0.0
    n_eff = 0.0
    """
    if how_many > 10000:
        galrange = [0,8885,how_many]
        steps = 2
        MAT = matrix(n.zeros(((7*(galrange[m+1]-galrange[m])),(7*(galrange[m+1]-galrange[m])))))
    if how_many < 10000:
        galrange = [0,how_many,how_many]
        MAT = matrix(n.zeros(((7*how_many),(7*how_many))))
        steps = 1
    """
    #for i in range(0,how_many):
    for i in range(0,n.shape(matches)[0]):
        gal1_index = matches[i][0]
        gal2_index = matches[i][1]
        gal1_bounds = bounds[gal1_index]
        gal2_bounds = bounds[gal2_index]
        #some of the stamps might be on edge of chip so the new stamp should not be too big or small
        xmin = n.min((gal1_bounds.getXMin(),gal2_bounds.getXMin()))
        if xmin < 1.0:
            xmin = 1

        ymin = n.min((gal1_bounds.getYMin(),gal2_bounds.getYMin()))
        if ymin < 1.0:
            ymin = 1

        xmax = n.max((gal1_bounds.getXMax(),gal2_bounds.getXMax()))
        if xmax > size:
            xmax = int(size)

        ymax = n.max((gal1_bounds.getYMax(),gal2_bounds.getYMax()))
        if ymax > size:
            ymax = int(size)
            #creat the stamp with appropriate size
        stampbox=galsim.BoundsI(int(xmin),int(xmax),int(ymin),int(ymax))
        newstamp = galsim.ImageD(stampbox,init_value = 0.0)
            #handle the diagonal terms by only populating diagonal part of MAT
        if gal1_index == gal2_index:
            cube = list()
            scidata = hdulist[(gal1_index+1)].data
            cube.append(scidata[0].flatten())
            cube.append(scidata[3].flatten())
            cube.append(scidata[4].flatten())
            cube.append(scidata[5].flatten())
            cube.append(scidata[6].flatten())
            cube.append(scidata[7].flatten())
            cube.append(scidata[8].flatten())
        if gal1_index != gal2_index:
            #find the overlap of each stamp with the new stamp
            #shift the galaxy bounds to correspond to the location of the new stamp
            overlap_1 = newstamp.bounds & gal1_bounds
            overlap_2 = newstamp.bounds & gal2_bounds
            #now put the stamp data into the new stamp
            #the first array of hdulist is empty
            scidata1 = hdulist[(gal1_index + 1)].data
            scidata2 = hdulist[(gal2_index + 1)].data
            #create and write the stamp data to a cube
            cube = list()
            #flux
            newstamp[overlap_1] += scidata1[0]
            newstamp[overlap_2] += scidata2[0]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
            #xc
            newstamp[overlap_1] += scidata1[3]
            newstamp[overlap_2] += scidata2[3]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
            #yc
            newstamp[overlap_1] += scidata1[4]
            newstamp[overlap_2] += scidata2[4]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
            #hlrd
            newstamp[overlap_1] += scidata1[5]
            newstamp[overlap_2] += scidata2[5]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
            #hlrb
            newstamp[overlap_1] += scidata1[6]
            newstamp[overlap_2] += scidata2[6]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
            #g1
            newstamp[overlap_1] += scidata1[7]
            newstamp[overlap_2] += scidata2[7]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
            #g2
            newstamp[overlap_1] += scidata1[8]
            newstamp[overlap_2] += scidata2[8]
            cube.append(newstamp.array.flatten())
            newstamp *= 0.0
        #make a holder that is size of the new stamp
        holder = n.copy(newstamp)*0.0
        #need to calculate the parameter crosscorelation
        #the variance in each pixel is equal to the sky_level + pixel flux
        #handle the diagonal terms by only populating diagonal part of MAT
        if gal1_index == gal2_index:
            for k in range(7):
                for j in range(7):
                    MAT[(gal1_index*7+k),(gal1_index*7+j)] = n.sum((cube[k]*cube[j])/(skylevel+cube[0]))
                    #now is a good time to find the null parameters (i.e. hlr_d hlr_b)
                    if k == j and MAT[(gal1_index*7+k),(gal1_index*7+j)] == 0.0:
                        #1 means delete
                        deletelist.append((gal1_index*7+k))
                        tokeepornot.append(1)
                    if k == j and MAT[(gal1_index*7+k),(gal1_index*7+j)] != 0.0:
                        tokeepornot.append(0)
        #the matches of different galaxies only require off-diagonal terms to be calculated in MAT
        if gal1_index != gal2_index:
            for k in range(7):
                for j in range(7):
                    MAT[(gal1_index*7+k),(gal2_index*7+j)] = n.sum((cube[k]*cube[j])/(skylevel+cube[0]))
                    MAT[(gal2_index*7+j),(gal1_index*7+k)] = n.sum((cube[k]*cube[j])/(skylevel+cube[0]))
    print deletelist
    print tokeepornot
    hdulist.close()
    #Remove all the rows and columns where hlrd and hlrb are zero
    #use numpy delete 'n.delete(matrix,list of row or column index,axis=0 or 1)'
    #axis = 0 is the row axis = 1 is the column
    MAT = n.delete(MAT,deletelist,axis=0)
    MAT = n.delete(MAT,deletelist,axis=1)

    #Invert the Fisher matrix to get to ellipticity variances
    IMAT = MAT.I
    #now we have the full matrix, contract it by removing the nuance parameters
    MAT_e = matrix(n.zeros(((2*how_many),(2*how_many))))
    eselect = list()
    counter = -2
    for k in range(n.size(tokeepornot)/7):
        if n.sum(tokeepornot[7*k:7*k+7]) == 1:
            #shift up a row/column
            counter += 6
            eselect.append(counter)
        if n.sum(tokeepornot[7*k:7*k+7]) == 0:
            #don't shift
            counter += 7
            eselect.append(counter)
    for k in range(n.size(eselect)):
        for j in range(n.size(eselect)):
            MAT_e[(2*k),(j*2)] = IMAT[eselect[k],eselect[j]]/N_exp
            MAT_e[(2*k+1),(j*2+1)] = IMAT[(eselect[k]+1),(eselect[j]+1)]/N_exp
            MAT_e[(2*k+1),(j*2)] = IMAT[(eselect[k]+1),eselect[j]]/N_exp
            MAT_e[(2*k),(j*2+1)] = IMAT[eselect[k],(eselect[j]+1)]/N_exp
    
    #Add in the shape noise to the diagonal terms
    
    for k in range((2*how_many)):
        for j in range((2*how_many)):
            if k == j:
               MAT_e[k,j] += sigma_SN**2
    #now invert the variance matrix and sum over each 2x2 to get Ce_inv
    SUBM = MAT_e.I
    for k in range(0,(2*how_many),2):
        for j in range(0,(2*how_many),2):
            
            if k == j:
                quad2 = SUBM[k,j]
                quad4 = SUBM[k+1,j+1]
                variance_m = ((MAT_e[k,j]-sigma_SN**2) + (MAT_e[k+1,j+1]-sigma_SN**2))/2.
                n_eff += sigma_SN**2/(sigma_SN**2 + (variance_m/N_exp))
            if k != j:
                quad2 = MAT_e[i,j]
                quad4 = MAT_e[i+1,j+1]
                thisn = sigma_SN**2/(((quad2+quad4)/2.0))
                n_eff += thisn
    return n_eff