#!/usr/bin/env python
# coding: utf-8


import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl

import pandas as pd
from scipy.spatial import distance_matrix


#sys.path.append("C:\Users\Cobal\Desktop\ZERNIKE\Codice_and_surface\Codice")

import ZernikeFunc as ZF
import SurfaceFunc as SF

import time 

SAVE_FIG = 0
PLOT = 0



Npixel = 25   # the plane edge in pixels...

Daa    = 3.    # the maximum distance between antigen and antibody binding site atoms..
Rs     = 6.    # the radius of the sphere that includes the patch..
Dpp = 1.       # the distance between points of the same patch (needed to remove islands)

ZOrder = 20   # the Zernike expansion order..

pdb_file= "..\RRM2\cluster2.dms"
verso = float(1)
respath = "..\RRM2\cluster2\positive"
name_ = "Cluster 2 of RRM2"
Npoint = int(1)

#sys.argv[1] = pdb_file
#sys.argv[2] = verso
#sys.argv[3] = respath
#sys.argv[4] = name_
#sys.argv[5] = Npoint

#if(len(sys.argv) != 6):
#	print("Error! Please insert:")
#	print("python SurfaceScreening_general dmsfile  verso  resdir name Npoint")
#	print("dmsfile, the surface file;")
#	print("verso, 1 for 'up' orientation of the patch, -1 for 'down';")
#	print("resdir, the directory where to save the outputs;")
#	print("name, a string to specify the output (i.e. pdb code);")
#	print("Npoint, the point interval for the screening, i.e. 1 for all points, 10 for one point every 10 points, etc.")
#	exit()

#pdb_file = sys.argv[1]
#verso = float(sys.argv[2])
#respath = sys.argv[3]
#name_ = sys.argv[4]
#Npoint = int(sys.argv[5])



if(verso == 1):
	ver = "up"
elif(verso == -1):
	ver = "down"

# loading surface file..

if(1):    
        
    
        surf_name =  pdb_file
    
        print("Processing protein %s with verso: (%s)\n"%(name_,verso))
        
        try:
            surf_ = pd.read_csv(surf_name)
        except:
            print("Given protein is not present!")
            exit()
    
        
    
    
        lag = len(surf_["x"])
        print(lag)

        surf = np.zeros((lag, 6))
        surf[:,:] = surf_[["x", "y", "z", "Nx", "Ny", "Nz"]]
    
    
        surf_obj = SF.Surface(surf[:,:], patch_num = 0, r0 = Rs, theta_max = 45)
        
        
        ### check sampled patch of real antibody bs
        if(PLOT):
            res1, c = SF.ConcatenateFigPlots([surf[:,:3]])
            SF.Plot3DPoints(res1[:,0], res1[:,1], res1[:,2], c, 0.3)


            
        ltmp = np.shape(surf)[0]

        
        index_possible_area = np.arange(ltmp)[::Npoint]  #MM]
     
        ltmp = len(index_possible_area)
        
        zernike_sampling_inv_a = np.zeros((121, ltmp))
    
        for i in range(ltmp):
            sys.stderr.write("\r Processing %i out of %i point"%(i, ltmp))
            sys.stderr.flush()
            # finding antigen patch, plane and zernike descriptors..
        
            patch, mask = surf_obj.BuildPatch(point_pos=index_possible_area[i], Dmin=.5)
            surf_obj.real_br = mask
        
            rot_patch, rot_ag_patch_nv = surf_obj.PatchReorientNew(patch, verso)
    
            
            z = surf_obj.FindOrigin(rot_patch)
            plane, weigths, dist_plane, thetas = surf_obj.CreatePlane(patch=rot_patch, z_c=z , Np=Npixel)
            new_plane = plane.copy()
            new_plane___ = plane.copy()
            
            
            if(np.shape(rot_patch)[1] == 4):
            
                new_plane_re = surf_obj.FillTheGap_everywhere(plane_=np.real(plane))
                new_plane_im = surf_obj.FillTheGap_everywhere(plane_=np.imag(plane))
            
                new_plane_re_ =  surf_obj.EnlargePixels(new_plane_re)
                new_plane_im_ =  surf_obj.EnlargePixels(new_plane_im)
            
                new_plane_ = new_plane_re_ + 1j*new_plane_im_/np.max(np.abs(new_plane_im_))
            else:
                new_plane = surf_obj.FillTheGap_everywhere(plane_=plane)
                ## enlarging plane..
                new_plane_ =  surf_obj.EnlargePixels(new_plane)
            
            
            try:
                zernike_env.img  = new_plane_
            except:
                zernike_env = ZF.Zernike2d(new_plane_)
            #br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
            br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)
            
            zernike_sampling_inv_a[:,i] = np.absolute(br_coeff)
            
        res_inv_ = np.row_stack([index_possible_area, zernike_sampling_inv_a])
    
        np.savetxt("%s/Selected_points_1out%d_zernike_invariant_%s_verso_%s.dat"%(respath, Npoint,name_, ver),res_inv_, fmt="%.4e")      
        np.savetxt("%s/Selected_points_1out%d_indexes_%s_verso_%s.dat"%(respath,Npoint, name_, ver), index_possible_area)




