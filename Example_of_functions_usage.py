#!/usr/bin/env python



import os, sys
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
import pandas as pd
from mayavi import mlab
import ZernikeFunc as ZF
import SurfaceFunc as SF


## Parameters definition
Npixel = 25   # the plane edge in pixels...
Rs     = 6    # the radius of the sphere that includes the patch..
Dpp = 1       # the distance between points of the same patch (needed to remove islands)
ZOrder = 20   # the Zernike expansion order..


surf_name_a = "../RRM2/cluster1.dms"

surf_a_ = pd.read_csv(surf_name_a)
    
l_a = len(surf_a_["x"])
print("Npoints", l_a)

surf_a = np.zeros((l_a, 6))
surf_a[:,:] = surf_a_[["x", "y", "z", "Nx", "Ny", "Nz"]]

### To inizialize the Surface class:
surf_a_obj = SF.Surface(surf_a[:,:], patch_num = 0, r0 = Rs, theta_max = 45)

### To plot surface in 3D:
res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)
        
## To isolate one patch of the surface:
## Example, the patch around point 5000..
patch, mask = surf_a_obj.BuildPatch(point_pos=5000, Dmin=.5)

## To rotate a patch:
rot_patch, rot_patch_nv = surf_a_obj.PatchReorientNew(patch, 1)


############ JUST PLOTS ################### 

##plotting surface + patch
res1, c = SF.ConcatenateFigPlots([surf_a_obj.surface[:,:3],patch[:,:3]])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)


##plotting patch + rotated patch
tmp1 = patch[:,:3] - np.mean(patch[:,:3], axis=0)
tmp2 = rot_patch[:,:3]

res1, c = SF.ConcatenateFigPlots([tmp1,tmp2])
SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)

##plotting rotated patch with normal vectors
SF.Plot3DPointsAndVectors(rot_patch[:,0], rot_patch[:,1], rot_patch[:,2],rot_patch_nv[:,0], rot_patch_nv[:,1], rot_patch_nv[:,2])

#############################################
    

## To project the patch on the xy plane...

z = surf_a_obj.FindOrigin(rot_patch)
plane, weigths, dist_plane, thetas = surf_a_obj.CreatePlane(patch=rot_patch, z_c=z , Np=Npixel)
            
new_plane = surf_a_obj.FillTheGap_everywhere(plane_=plane)
## enlarging plane..
new_plane_ =  surf_a_obj.EnlargePixels(new_plane)


# plotting res..
fig, ax = mpl.subplots(1,2, dpi = 150)
ax[0].imshow(plane)
ax[1].imshow(new_plane_)
ax[0].axis('off')
ax[1].axis('off')
ax[0].set_title("original")
ax[1].set_title("processed")
mpl.show()


# # Zernike representation

### To initialize the Zernike2D class:

zernike_env = ZF.Zernike2d(new_plane_)



### To compute and visualize a complex Zernike polynom...

n = 6
m = 2

Z_nm_1 = zernike_env.ComputeMoment(n,m)

fig, ax = mpl.subplots(1,2, dpi = 150)
ax[0].imshow(Z_nm_1.real)
ax[1].imshow(Z_nm_1.imag)
ax[0].axis('off')
ax[1].axis('off')
ax[0].set_title("Re[Z(%d,%d)]"%(n,m))
ax[1].set_title("Im[Z(%d,%d)]"%(n,m))
mpl.show()


#zernike_env.ComputeCoeff_nm(zernike_env.img, 2,0)



#### To decompose an image in the Zernike basis:
recon, coeff = zernike_env.ZernikeRecostruction(order=ZOrder, PLOT=0)
   
### ZernikeRecostruction returns also the reconstructed plane, so it is slow. Use ZernikeDecomposition, which returns only coeff but it is much faster.



        

#### To obtain the invariant descriptors:
## Coefficients are complex numbers.
## Taking the modulus, you obtain real numbers, which are rotationally invariant..

coeff_inv = np.absolute(coeff)




# plotting res..
fig, ax = mpl.subplots(1,3, dpi = 150)
ax[0].imshow(plane)
ax[1].imshow(new_plane_)
ax[2].imshow(recon.real)
ax[0].axis('off')
ax[1].axis('off')
ax[2].axis('off')
ax[0].set_title("original")
ax[1].set_title("processed")
ax[2].set_title("reconstructed")



mpl.figure(dpi = 150)
mpl.plot(coeff_inv)
mpl.xlabel("index")
mpl.ylabel("Invariant")
mpl.show()