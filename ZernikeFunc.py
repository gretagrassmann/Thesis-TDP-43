#!python

import os, sys
import numpy as np
import matplotlib.pyplot as mpl
import scipy as sp
import scipy.ndimage as spd
from matplotlib.image import imread
from matplotlib import cm



try:
    from mayavi import mlab
    MLAB_LAB = 1
  
except:
    print("mayavi not present")
    MLAB_LAB = 0
    

from scipy.special import sph_harm

def myflip(m, axis):
    if not hasattr(m, 'ndim'):
        m = asarray(m)
    indexer = [slice(None)] * m.ndim
    try:
        indexer[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError("axis=%i is invalid for the %i-dimensional input array"
                         % (axis, m.ndim))
    return m[tuple(indexer)]


def IsolateIsosurface(myprot, minV, maxV):
    #### This function groups points nearer than minD..
    
    DEB_ = 0
    
    lx, ly, lz = np.shape(myprot)
    
    prot = np.copy(myprot)
    mask = np.logical_and(prot>=minV, prot<= maxV)
    #print(mask)
    prot[:,:,:] = 0
    prot[mask] = 1.
    prot_label = np.copy(prot)
    #print(np.unique(prot, return_counts=True))
    #print(len(mask))
    
    
    # starting from label = 2
    lab = 2
    
    #surf_label[0] = lab
    #surf_tmp[0] = lab
    
    ## defining probe points..
    tmp = np.zeros((3,3,3))
    x, y, z = np.where(tmp==0)
    x = x-1
    y = y-1
    z = z-1
    
    #pos_s = np.where(surf_tmp == lab)
    
    # computing number of points without label..
    Nleft = np.sum(prot != 0)

    #print(Nleft)
    # starting iterating over different surfaces..
    while(Nleft > 0):
        
        count = 1
        
        pos__ = np.where(prot != 0)
        if(DEB_):
            print("pos__",np.shape(pos__))
        
        # seeding: first unlabeled point takes lab label.. 
        #print("---",prot_label[pos__[0][0]])
        prot_label[pos__[0][0], pos__[1][0], pos__[2][0]] = lab
        prot[pos__[0][0], pos__[1][0], pos__[2][0]] = lab
 
        if(DEB_):
            print("prot", prot)
        

        # iterating to find points belonging to the same surface...
        while(count > 0):
            
            #print(count)
            count = 0
            pos_s = np.where(prot == lab)
            #print("lll",len(pos_s))
            if(len(pos_s) == 0):
                break
            # creating mask for points still to be processed..
            mask = np.logical_and(prot > 0, prot != lab)
            if(DEB_):
                print("l", np.shape(pos_s))
            for i in range(np.shape(pos_s)[1]):
                if(DEB_):
                    print("pos",i,np.shape(pos_s) )
                
                if(np.shape(pos_s)[1] == 1):
                    xxxx =  pos_s[0]
                    yyyy =  pos_s[1]
                    zzzz =  pos_s[2]
                else:
                    
                    xxxx = pos_s[0][i]
                    yyyy = pos_s[1][i]
                    zzzz = pos_s[2][i]
                x_ = x + xxxx
                y_ = y + yyyy
                z_ = z + zzzz
                
                mask_x = np.logical_and(x_ >= 0, x_<lx )
                mask_y = np.logical_and(y_ >= 0, y_<ly )
                mask_z = np.logical_and(z_ >= 0, z_<lz )
                
                mask_xyz = np.logical_and(mask_x, np.logical_and(mask_y, mask_z))
                
                x_ = x_[mask_xyz]
                y_ = y_[mask_xyz]
                z_ = z_[mask_xyz]
                if(DEB_):
                    print("xyz",x_, y_,z_)
                
                for j in range(len(x_)):
                    if(prot[x_[j], y_[j], z_[j]] == 1):
                        prot[x_[j], y_[j], z_[j]] = lab
                        prot_label[x_[j], y_[j], z_[j]] = lab
                        count += 1

                if(DEB_):
                    print("prot_2", prot)

                # removing processed point from system..
                prot[xxxx,yyyy,zzzz] = 0
                mm = np.logical_and(prot > 0, prot != lab)
       
                #print("00", np.sum(surf_tmp == 0),"out of", len(surf_tmp), "lab", lab)
                
                # creating mask for  point still to be processed..
                mmm = np.logical_and(mm,mask)
            
                mask = np.logical_and(prot > 0, prot != lab)
                if(DEB_):
                    print("prot_c", prot)

        
        # creating a new label..
        lab += 1
        # looking for how many points still to be processed...
        Nleft = np.sum(prot != 0)
        sys.stderr.write("\rleft %d"%Nleft)
        sys.stderr.flush()
    return(prot_label)




def log10_factorial(n):
    '''
    This function recursively computes  the log10 of a factorial.
    
    Input:
    - n, the interger to compute n!
    
    Output:
    - log10(n!)
    '''
    
    if(n <=1):
        return(0)
    else:
        return(np.log10(n)+ log10_factorial(n-1))




class Zernike2d:
    '''
    This class performs the 2D decomposition of a figure in its Zernike descriptors
    '''



    def __init__(self, imagefile):

        if(type(imagefile) == str):

            self.img = self.PrepareImage(imagefile)
        else:
            self.img = imagefile


        Nl  = np.shape(self.img)[0]
        self.Nl = Nl
        #come_back = np.zeros((Nl,Nl), dtype=complex)

        self.x, self.y = self.BuildPlane(Nl)
        self.r, self.t = self.FromCartesianToPolarPlane(self.x,self.y)

        tmp = np.ones(np.shape(self.img))
        tmp = self.CircleImage(tmp)
        self.npix = np.sum(tmp)

        self.zernike_dict = {}



    def CircleImage(self,image):

        l,tmp = np.shape(image)

        new_image = image.copy()

        r  = ((l-1)/2.)
        r2 = r**2
        origin = np.array([r+1, r+1])


        for i in range(l):
            for j in range(i,l):
                d2 = (i-r)**2 + (j-r)**2
                if(d2>r2):
                    new_image[i,j] = 0
                    new_image[j,i] = 0
        return(new_image)

    def PrepareImage(self, datafile):
        data = imread(datafile)
        data = data[:,:,0]
        CUT = 1
        x,y = np.shape(data)
        l = np.min([x,y])
        if(l%2 == 0):
            l -= 1
        new_image = np.zeros((l,l))

        r  = ((l-1)/2.)
        r2 = r**2
        origin = np.array([r+1, r+1])

        if(x < y):

            start = int((y-l)/2.)
            new_image[:,:] = data[:l,start:start+l]
        elif(y < x):

            start = int((x-l)/2.)
            new_image[:,:] = data[start:start+l,:l]
        else:
            new_image[:,:] = data[:l,:l]

        if(CUT):
            for i in range(l):
                for j in range(i,l):
                    d2 = (i-r)**2 + (j-r)**2
                    if(d2>r2):
                        new_image[i,j] = 0
                        new_image[j,i] = 0
        return(new_image)

    def ComputeDot(self,A, B):



        c= np.sum(A*np.conjugate(B))/float(self.npix)
        return(c)

    def ComputeCoeff_nm(self,F, n,m):
        Nl, tmp = np.shape(F)
        dx = 1./(Nl-1)

        #x, y = self.BuildPlane(Nl)
        #r, t = self.FromCartesianToPolarPlane(x,y)
        Z = self.ComputeMoment(n,m)

        c = self.ComputeDot(F,Z)*float(n+1)
        return(c)

    def FromPolarToCartesian(self,r,theta):
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        return(x,y)

    def FromCartesianToPolarPlane(self,x,y):

        l, tmp = np.shape(x)
        R_ = np.zeros((l,l))
        Theta_ = np.zeros((l,l))

        for i in range(l):
            for j in range(l):

                r, t = self.FromCartesianToPolar(x[i,j],y[i,j])
                R_[i,j] = r
                Theta_[i,j] = t
        return(R_, Theta_)

    def FromCartesianToPolar(self, x,y):

        r = np.sqrt(x**2 + y**2)

        if(y==0 and x >0):
            theta = 0
        elif(y==0 and x <0):
            theta = np.pi
        else:
            t = np.arctan(np.abs(y/x))

            if(x> 0 and y>0):
                theta = t
            elif(x<0 and y<0):
                theta = t + np.pi
            elif(x<0 and y>0):
                theta = np.pi - t
            elif(x>0 and y<0):
                theta = 2*np.pi - t
            elif(x==0 and y >0):
                theta = np.pi/2.
            elif(x==0 and y <0):
                theta = 3.*np.pi/2.
            else:
                theta = 0.
        return(r, theta)


    def R_nm(self,n,m, Lr):
        '''
        This function computes the radial (R_nm) term 
        of the Zernike polynomial.
        Input:
        n, the order;
        m, the angular coefficient;
        Lr, the discretization lenght;

        Output:
        R_nm_, the radial matrix in the plane [-1,1]
        '''
        
        rr = self.r.copy()
        mask = rr == 0
        rr[mask] = 1

        log10_r = np.log10(rr)

        R_nm_ = np.zeros(Lr)
        R_nm_0 = 0
        if((n-m)%2 != 0):
            return(R_nm_)
        else:
            diff= int((n-m)/2.)
            summ = int((n+m)/2.)
            for l in np.arange(0, diff+1):
                ## using log for product..
                num = log10_factorial(n-l) + (n-2.*l)*log10_r
                den = log10_factorial(l) + log10_factorial(summ - l)  + log10_factorial(diff - l)

                if(n-2.*l == 0):
                    num0 = log10_factorial(n-l)
                    R_nm_0 += (-1)**l*10.**(num0-den)
                
                R_nm_ += (-1.)**l*10.**(num - den)
            R_nm_[mask] = R_nm_0
            return(R_nm_)

    def Phi_m(self, theta, m):
        '''
        This function computes the complex e^(i*m*theta) term 
        of the Zernike polynomial.
        '''
        phi = np.cos(m*theta) + 1j*np.sin(m*theta)
        return(phi)



    def BuildPlane(self,N):
        '''
        This function intitalizes a plane in [-1,1], with a discretization of NxN points. 
        '''
        
        plane_x = np.zeros((N,N))
        plane_y = np.zeros((N,N))

        Nr = int((N-1)/2.)
        dx = 1./Nr
        x = np.arange(0,N)*dx - 1.
        x_f = myflip(x,axis=0)
        for i in range(N):
            plane_x[i,:] = x
            plane_y[:,i] = x_f
        return(plane_x, plane_y)



    def CountMoment(self, n):
        '''
        This function computes the number of moment that an expansion to the n order will produce.
        '''

        if(n%2 == 0):
            N_z = n/2 + 1
            for i in range(1, int(n/2.)+1):
                N_z += 2*i
        else:
            N_z = 0
            for i in range(1, int((n+1)/2.)+1):
                N_z += 2*i
        return(N_z)



    def ComputeMoment(self,n,m):
        '''
        This function computes the Z_nm moment.
        Input:
        n, the order parameter;
        m, the angular parameter (<=n).

        Output:
        Z_nm, a complex matrix of the polynomial.
        '''
        
        Z_nm_ = self.zernike_dict.get((n,m))

        if Z_nm_ is None:
            Z_nm = np.zeros((self.Nl,self.Nl), dtype=complex)

            r_part = self.R_nm(n,np.abs(m),[self.Nl,self.Nl])

            Z_nm = r_part*self.Phi_m(self.t, m)

            Z_nm[np.isnan(Z_nm)] =0
            Z_nm_ = self.CircleImage(Z_nm)  #normalization factor.. #*np.sqrt((n+1))
            #Z_nm_ /= np.sqrt(np.sum(Z_nm_*Z_nm_*dx**2))

            self.zernike_dict[(n,m)] = Z_nm_
            return(Z_nm_)
        else:
            return(Z_nm_)





    def ZernikeRecostruction(self, order, PLOT = 1):
        '''
        This function performes the decomposition and recostruction of the function (image) in the Zernike basis. 
        USE "ZernikeRecostruction", which is FASTER if only interested in decomposition.
        
        Input:
        -order, the maximum order of the expansion (N)
        -PLOT, 0 (1) to plot (not plot) the result.
 
        Output:
        - come_back, a matrix with the recomposed image. 
        - c_list, a list of complex zernike coefficients (c_nm).
        
        Notes:
        1) The decomposed image is self.img
        2) The order of the coefficients is: for increasing n, alle possible m in increasing order). Es. (n=0,m=0), (n=1, m=1)...
        '''

        ROT_A = 0, #np.pi*30/180.
        c_list = []
        come_back = np.zeros((self.Nl, self.Nl), dtype=complex)
     
        ## cycling over index n (order)..
        for n in range(order+1):
            ### cycling over moment m in (0, n)
            for m in range(0, n+1):
                if((n-m)%2 == 0):
                    c = self.ComputeCoeff_nm(self.img, n,m)
                    c_list.append(c)
         
                    moment_ = self.ComputeMoment(n,m)

                    c_abs = np.absolute(c)
                    c_phi = np.arctan2(c.imag, c.real)

                    m_abs = np.absolute(moment_)
                    m_phi = np.arctan2(moment_.imag, moment_.real)


                    #come_back += c*moment_
                    if(m != 0):
                        come_back += c_abs*m_abs*np.exp(1j*(m*(m_phi/float(m) + ROT_A) + c_phi))
                    else:
                        come_back += c_abs*m_abs*np.exp(1j*(m_phi + c_phi))
     
                    if(PLOT):
                        fig, ax = mpl.subplots(1,2, dpi = 150)
                        ax[0].imshow(moment_.real)
                        ax[1].imshow(moment_.imag)
                        mpl.show()

        if(PLOT):
            fig, ax = mpl.subplots(1,2, dpi = 150)
            ax[0].imshow(come_back.real)
            ax[1].imshow(come_back.imag)
            mpl.show()

        return(come_back, c_list)


    def ZernikeDecomposition(self, order):
        '''
        This function performes the decomposition of the function (image) in the Zernike basis.
        
        Input:
        -order, the maximum order of the expansion (N)
        
        Output:
        - c_list, a list of complex zernike coefficients (c_nm).
        
        Notes:
        1) The decomposed image is self.img
        2) The order of the coefficients is: for increasing n, alle possible m in increasing order). Es. (n=0,m=0), (n=1, m=1)...
        '''
    
        c_list = []
     
        ## cycling over index n (order)..
        for n in range(order+1):
            ### cycling over moment m in (0, n)
            for m in range(0, n+1):
                if((n-m)%2 == 0):
                    c = self.ComputeCoeff_nm(self.img, n,m)
                    c_list.append(c)
        
        return(c_list)

