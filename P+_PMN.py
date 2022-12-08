import matplotlib as plt
import numpy as np
import astropy
from astropy import units as u

import astropy.coordinates as coord
from astropy.io import ascii
import astropy.constants as const
from astropy.coordinates import SkyCoord
import emcee
import corner
    
from astropy.cosmology import wCDM
from astropy.cosmology import w0wzCDM
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import w0waCDM 

from time import time
# from Cosmography_Class import Cosmography
import sys
import pymultinest as pmn





###########
#LOAD IN PANTHEON+
direc = 'DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/'#change to your own directory of data
table = ascii.read(direc+'Pantheon+SH0ES.dat')
tot_SNe= len(table)
#################
#run with stat only or with systematic covariance matrix

model = sys.argv[1]
# Cosmography_QUAD_DIP'
z_frame = sys.argv[2]
# z_frame = 'zCMB'
z_low = 0.023
z_up = 0.8
nlp=4000 #number of live points to sample over
x = 0.1 #decay scale
S_q = x/np.log(2)

#stuff for saving chain files
chains_str = "PMN_set_S_q/"+model+"_"+z_frame+"nlp="+str(nlp)+ '_'+str(x)+'np.log(2)_'+str(z_low)+'<z'
outputfiles_basename= chains_str #directory to save chains

withsys = True
withcal= False

if withsys:
    print("With Systematics")
    statsys = np.loadtxt(direc+'Pantheon+SH0ES_STAT+SYS.cov', skiprows =1)
    cov = np.reshape(statsys, (tot_SNe, tot_SNe))
    chains_str += "WithSys_"
    
else:
    print('Without systmatics')
    stat = np.loadtxt(direc+'Pantheon+SH0ES_STATONLY.cov', skiprows=1)
    tot_SNe = len(table)
    cov = np.reshape(stat, (tot_SNe, tot_SNe))
    chains_str += "WithCal_"
    
if withcal:
    print('With cal')
    print('calibrator not yet coded')
    data = table
    C = cov
    raise NameError
else:
    print('Without calibrators')
    data = table[(table['IS_CALIBRATOR']==0)& (table['zHD']>= z_low) &(table['zHD']<=z_up) ]
    indices = np.where((table['IS_CALIBRATOR']==0)&(table['zHD']>= z_low)&(table['zHD']<= z_up))[0]#finding indices for non-calibrators 
    C = cov[indices][:,indices] #using covariance matrix without calibrators
    # print(C.shape)


#load in coordinate system

ra = (data['RA']); dec = data['DEC']

#load in redshifts and needed columns
if z_frame == 'zHEL':
    z = data['zHEL']
elif z_frame == 'zCMB':
    z = data['zCMB']
elif z_frame == 'zHD':
    z = data['zHD']    
else:
    print('specify redshift frame!')
    raise NameError

MU_SH0ES = data['MU_SH0ES']; MU_SH0ES_ERR_DIAG = data['MU_SH0ES_ERR_DIAG']
m_b_corr = data['m_b_corr']; m_b_corr_err_DIAG = data['m_b_corr_err_DIAG']


class Cosmography:
    '''Class for calculating distance luminosity relations'''

    def __init__(self, zarr, m_b, covmat, ra, dec):
     
        self.covmat = covmat
        self.invC = np.linalg.inv(covmat)
        self.z = zarr
        self.m_b = m_b
        #CMB  DIPOLE COORDINATES
        self.CMB = SkyCoord(l=264.021, b=48.253, unit = 'degree',frame='galactic')
        self.c = const.c* u.s/u.km 
        self.Sne = SkyCoord(ra, dec, unit='degree', frame='icrs')
        #coordinates of quadrupole from ***
        self.n1 = SkyCoord(l=118, b=85, frame = 'galactic', unit='deg')
        self.n2 = SkyCoord(l = 341, b =4, frame = 'galactic', unit='deg')
        self.n3 =SkyCoord( l = 71, b = - 4, frame = 'galactic', unit='deg')

        self.theta1 = self.n1.separation(self.Sne).radian
        self.theta2= self.n2.separation(self.Sne).radian
        self.theta3 = self.n3.separation(self.Sne).radian
        self.e_dip = self.CMB.separation(self.Sne).radian
        print(self.theta1)


    def F_decay(self, z, S):
        '''exponential decay function for Q dipole'''
        F = np.exp(-z/S)
        return(F)

    # def sep(self, lat2, lon2):
    #     '''calculates angular separation between two points'''
    #     lat1 = self.dec_SNe; lon1 = self.ra_SNe;
    #     sdlon = np.sin(lon2 - lon1)
    #     cdlon = np.cos(lon2 - lon1)
    #     slat1 = np.sin(lat1)
    #     slat2 = np.sin(lat2)
    #     clat1 = np.cos(lat1)
    #     clat2 = np.cos(lat2)

    #     num1 = clat2 * sdlon
    #     num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    #     denominator = slat1 * slat2 + clat1 * clat2 * cdlon
    #     e =  np.arctan2(np.hypot(num1, num2), denominator)
    #     return(e)
    # def sep(self, lat2, lon2):
    #     '''calculates angular separation between two points'''
    #     lon1 = self.dec_SNe; lat1 = self.ra_SNe;
    #     sdlon = np.sin(lon2 - lon1)
    #     cdlon = np.cos(lon2 - lon1)
    #     slat1 = np.sin(lat1)
    #     slat2 = np.sin(lat2)
    #
       
    def prior(self, cube, npar, ncube):
        if model == 'DIP_QUAD':
        #h0_mono, M, q0, j, L1, L2, S= theta #j=Omega_k0 - j0

            cube[0] = cube[0] * 50. + 50. #H0
            cube[1] = cube[1] * 30. - 20. #M
            cube[2] = cube[2] * 20. - 10. #qm
            cube[3] = cube[3] * 30. - 10.#qd
            cube[4] = cube[4] * 20 - 10. #j
            cube[5] = cube[5] * 4. -2. #l1
            cube[6] = cube[6] * 4. -2. #l2
            cube[7] = cube[7] * 2.0 +0.01 #s1
            # cube[8] = cube[8] * 2. +0.01 #S2
        elif model == 'DIP_ONLY':
            # H0, M, qm, qd, j, S = theta #j=Omega_k0 - j0
            cube[0] = cube[0] * 50. + 50. #H0
            cube[1] = cube[1] * 30. - 20. #M
            cube[2] = cube[2] * 20. - 10. #qm
            cube[3] = cube[3] * 30. - 10.#qd
            cube[4] = cube[4] * 20 - 10. #j
            cube[5] = cube[5]* 2.0 +0.01 #sdip
        elif model == 'QUAD_ONLY':
            # H0_mono, M, q0, j, L1, L2, S= theta #j=Omega_k0 - j0
            cube[0] = cube[0] * 50. + 50.
            cube[1] = cube[1] * 30. - 20.
            cube[2] = cube[2] * 20. - 10. #q0
            cube[3] = cube[3] * 20. - 10.#J
            cube[4] = cube[4] * 4. -2. #L1
            cube[5] = cube[5] * 4. -2. #l2
            # cube[6] = cube[6] *2 + 0.001 #S_Q
        else:
            raise NameError
    def llhood(self, model_param, ndim, ncube):
        #fitting only for quadrupole in the hubble parameter
        z = self.z
        m_b = self.m_b
        if model== 'DIP_QUAD':
            H0_mono, M, qm, qd, j, L1, L2, S1, = [model_param[i] for i in range(8)] #j=Omega_k0 - j0
            theta1 = self.theta1
            theta2 = self.theta2
            theta3 = self.theta3
            e = np.cos(self.e_dip)
            # theta1, theta2, theta3 = self.quadrupole_coord()
            S2 = S_q
            q0 = qm + e*qd *self.F_decay(z, S1) #equation of dipolex
            h0_aniso = L1*(np.cos(theta1))**2 + L2*(np.cos(theta2))**2 -(L1 + L2)*(np.cos(theta3))**2 
            H0 = H0_mono*(1+h0_aniso*(self.F_decay(z,S2)) )
        elif model=='QUAD_ONLY':
            H0_mono, M, q0, j, L1, L2,   = [model_param[i] for i in range(6)]
            # theta2, theta3 = self.quadrupole_coord()
            theta1 = self.theta1
            theta2 = self.theta2
            theta3 = self.theta3
            # H = 1/3 * thet - e_mu*e_nu*sigma
            S2 = S_q
            # print(theta1, theta2, theta3)
            # print(np.cos(theta1)**2)
            H0 = H0_mono*(1 +(L1*(np.cos(theta1))**2 + L2*(np.cos(theta2))**2 -(L1 + L2)*(np.cos(theta3))**2 )* self.F_decay(z,S2)) #dipole equation
            
        elif model =='DIP_ONLY':
             #anisotropy parameters: qm: monopole value, qd: dipole value, S: dipole decay scale
            H0, M, qm, qd, j, S=  [model_param[i] for i in range(6)]  #j=Omega_k0 - j0
            # e = np.cos(self.CMB.separation(self.Sne).radian) #dot product
            e = np.cos(self.e_dip)
            # e = np.cos(astropy.coordinates.angular_separation(self.ra_CMB,self.dec_CMB, self.ra_SNe, self.dec_SNe))
            # e = np.cos(self.sep(self.dec_CMB, self.ra_CMB))
            q0 = qm + e*qd *self.F_decay(z, S) #equation of dipole
        else:
            raise NameError

        

        dL1 = 1/H0
        dL2 = (1-q0)/(2*H0)
        dL3 = (-1 + 3*(q0**2) +q0 + j)/(6*H0)
        dL = dL1*z + dL2*(z**2) + dL3*(z**3) #equation 2 in hayley paper
        mod = 5*np.log10(self.c*dL) + 25 + M
        delta_m =m_b - mod
        chi2 = np.matmul(np.transpose(delta_m), np.matmul(self.invC, delta_m)) #equation8
        return -0.5*chi2
    def run(self, nlp):
        print('running chains!')
        '''runs MCMC for given model with N_steps, takes z array, array of mb'''
        if model== 'DIP_QUAD':
            npar = 8
        elif model== 'QUAD_ONLY':
            npar = 6
        elif model== 'DIP_ONLY':
            npar = 6
        else:
            raise NameError
        pmn.run(self.llhood, self.prior, npar, verbose=False, n_live_points=nlp, outputfiles_basename=chains_str+'-')
st = time()



cosmo = Cosmography(z, MU_SH0ES, C, ra, dec)
cosmo.run(nlp)
end = time()
dur = (end - st) / 60.
print("The total time it took is ", dur, " minutes")
