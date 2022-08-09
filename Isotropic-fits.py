import matplotlib as plt
import numpy as np
import astropy
from astropy import units as u
from astropy.units import cds
cds.enable()
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import ascii
import emcee
import corner
    
from astropy.cosmology import wCDM
from astropy.cosmology import w0wzCDM
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import w0waCDM 

direc = '/data/jac304/DataRelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/'#change to your own directory of data

table = ascii.read(direc+'Pantheon+SH0ES.dat')
stat = np.loadtxt(direc+'Pantheon+SH0ES_STATONLY.cov', skiprows=1)
dstat = np.reshape(stat, (1701, 1701))
statsys = np.loadtxt(direc+'Pantheon+SH0ES_STAT+SYS.cov', skiprows =1)
csys = np.reshape(statsys, (1701, 1701))

data = table[table['IS_CALIBRATOR']==0] #selecting without calibrators for now

indices = np.where(table['IS_CALIBRATOR']==0)[0]#finding indices for non-calibrators 
print(dstat.shape)
C= dstat[indices][:,indices] #using covariance matrix without calibrators
print(C.shape)
print(C)
invC = np.linalg.inv(C) #inverse covariance matrix
print(len(C))

#load in coordinate system
ra = (data['RA']); dec = data['DEC']
c = astropy.coordinates.SkyCoord(ra = ra*u.degree, dec = dec*u.degree)

#load in redshifts and needed columns
zhel = data['zHEL']
zCMB = data['zCMB']; zCMBERR = data['zCMBERR']

MU_SH0ES = data['MU_SH0ES']; MU_SH0ES_ERR_DIAG = data['MU_SH0ES_ERR_DIAG']
m_b_corr = data['m_b_corr']; m_b_corr_err_DIAG = data['m_b_corr_err_DIAG']




galactic = c.galactic #use galactic coordinates 
gal = SkyCoord(ra[:], dec[:], frame='fk5', unit=u.deg)

def integrand(z, M, H0, Om0):
    '''integrand part of luminoscity distance relation to pass to the scipy integrate.quad function'''
    omega_l = 1. - Om0
    omega_k = 0.
    w = -1
    E = np.sqrt(Om0*(1+z)**3+ omega_l*(1+z)**(3*(1+w))+(omega_k)*(1+z)**2)
    return (1/E)

def integral(z_arr,  M, H0, Om0):
    '''integral part of distance luminoscity function'''
   #parameters to sample over
    args = ( M, H0, Om0 )
    omega_l = 1. - Om0
    I_arr = np.empty(shape=z_arr.shape)
    for z in z_arr:
        i = np.where(z_arr == z)
        I = scipy.integrate.quad(integrand,0, z, args = args)[0] 
        I_arr[i] = I
    return(I_arr)

def log_likelihood(theta, z, m_b_corr,model):
    '''log likelihood function m '''
    if model =='FlatLambdaCDM':
        M, H0, Om0 = theta
        w = -1
    #     dL = c/H0 *(1+z)*integral(z, M, Om0, H0)
    #     omega_l = 1 - omega_m #assuming flat
        flc = FlatLambdaCDM(H0, Om0)
        mod = (flc.distmod(z).value) + M 

    elif model == 'wCDM' :
        M, H0, Om0, w0 = theta
        Ode0 = 1 - Om0
        flc = wCDM(H0, Om0, Ode0, w0)
        mod = flc.distmod(z).value + M 
        
    elif model =='Cosmography':
#         H0, M, qm, qd, j, S1, L1, L2, S2= theta #j=Omega_k0 - j0
#         H0, M, qm, qd, j, S= theta #j=Omega_k0 - j0
#####################
# dipole only in q
#         e = sep(dec_CMB,ra_CMB)
    
#         q0 = qm + e*qd *F_decay(z, S)
    
#         q0 = qm + e*qd *F_decay(z, S1)
#         H = 1/3*thet - e_mu*e_nu*sigma
#only quadrupole in H

#         H0_mono, M, q0, j, L1, L2, S2= theta #j=Omega_k0 - j0
#         H0 = H0_mono*(1 +(L1*(np.cos(theta1))**2 + L2*(np.cos(theta2))**2 -(L1 + L2)*(np.cos(theta3))**2 )* F_decay(z,S2))
 

        
# dipole in q & quadrupole in H
        H0_m, M, qm, qd, j, S1, L1, L2, S2= theta #j=Omega_k0 - j0
        q0 = qm + e*qd *F_decay(z, S1)
   
        H0 = H0_m*(1 +(L1*(np.cos(theta1))**2 + L2*(np.cos(theta2))**2 -(L1 + L2)*(np.cos(theta3))**2 )* F_decay(z,S2))
 #SAME FOR ALL
        dL1 = 1/H0
        dL2 = (1-q0)/(2*H0)
        dL3 = (-1 + 3*(q0**2) +q0 + j )/(6*H0)
        dL = dL1*z + dL2*(z**2) + dL3*(z**3) #equation 2 in hayley paper
#             dL[i] = dl
        mod = 5*np.log10(c*dL) + 25 + M
    elif model =='FLRW_Cosmography':
        
        H0, M, q0, j = theta #j=Omega_k0 - j0
#         dL = np.empty(shape=z.shape)
#         z_arr = z
#         for z in z_arr:
#             i = np.where(z_arr == z)
#             dL1 = 1/H0
#             dL2 = (1-q0)/ (2*H0)
#             dL3 = (-1 + 3*q0**2 +q0 + j )/(6*H0)
#             dl = dL1*z + dL2*z**2 + dL3*z**3 #equation 2 in hayley paper
#             dL[i] = dl
#         mod = 5*np.log10(dL) + 25 + M

        dL1 = 1/H0
        dL2 = (1-q0)/(2*H0)
        dL3 = (-1 + 3*(q0**2) +q0 + j )/(6*H0)
        dL = dL1*z + dL2*(z**2) + dL3*(z**3) #equation 2 in hayley paper
        mod = 5*np.log10(c*dL) + 25 + M
# #          mod = (dL)+ M
#         dL1 = 1/H0
#         dL2 = (1-q0)/ (2*H0)
#         dL3 = (-1 + 3*q0**2 +q0 + j )/(6*H0)
#         dl = dL1*z + dL2*z**2 + dL3*z**3 
#         mod = 5*np.log10(dL) + 25 + M

  
    elif model =='w0waCDM':
        M, H0, Om0, w0, wa = theta
        Ode0 = 1 - Om0
        flc = w0waCDM(H0, Om0, Ode0, w0, wa)
        mod = flc.distmod(z).value + M 
    
  
    elif model =='w0wzCDM':
        M, H0, Om0, wz = theta
        Ode0 = 1 - Om0
        flc = w0wzCDM(H0, Om0, Ode0, w0, wz)
        mod = flc.distmod(z).value + M 
   
    else:
        raise NameError
        print('Invalid model, valid models are: wCDM, w0waCDM, cosmography , FlatLambdaCDM')
    
    
    delta_m =m_b_corr - mod
    chi2 = np.matmul(np.transpose(delta_m), np.matmul(invC, delta_m)) #equation8
    L = -0.5*chi2
    
    if np.isnan(L):
        return -  np.inf
    
    return(L)

def log_prior(theta, model):
    '''priors on sampled parameters for each model, returns -inf if outside of priors, 0 if inside'''
    if model =='FlatLambdaCDM':
        M, H0, Om0 = theta
        Ode0 = 1-Om0
        w0 = -1 
        if 0.0 <= Om0 < 1.0  and -20. < M  < 20. and  20. < H0 < 100.:
            return 0.0
        else:
            return -np.inf
        
    elif model =='wCDM':
        M, H0, Om0, w0, = theta
        if 0.0 <= Om0 < 1.0  and -20. < M  < 20. and  20. < H0 < 100. and -5 < w0 < 0.:
            return 0.0
        else :
            return -np.inf
    elif model =='w0waCDM':
        M, H0, Om0, w0, wa = theta
        if 0.0 <= Om0 < 1.0  and -20. < M  < 20. and  20. < H0 < 100. and -5 < w0 < 0.:
            return 0.0
        else:
            return -np.inf
    elif model =='w0wzCDM':
        M, H0, Om0, w0, wa = theta
        if 0.0 <= Om0 < 1.0  and -20. < M  < 20. and  20. < H0 < 100. and -5 < w0 < 0.:
            return 0.0
        else:
            return -np.inf
        

#     elif model =='FLRW_Cosmography':
#         H0, M, q0, j = theta #j=Omega_k0 - j0
#         if  -20. < M  < 10. and  20. < H0 < 100. and -5. < j < 5. and -5. < q0 < 5.:
#             return 0.0    
#         else:
#             return -np.inf
    
#     elif model =='Cosmography':
        
#         #dipole in q
# #         H0, M, qm, qd, j, S = theta #j=Omega_k0 - j0

#         #quadrupole in H
# #         H0_mono, M, q0, j, L1, L2, S= theta #j=Omega_k0 - j0
#         #quadrupole in H and dipole in q 
#         H0_m, M, qm, qd, j, S1, L1, L2, S2 = theta #j=Omega_k0 - j0
#         if  -20. < M  < 10. and  20. < H0_m < 100. and -5. < j < 5  and  0.01 < S1 < 0.1 and  0.01 < S2 < 0.1 and -2. < L1 < 2. and -2. < L2 < 2.:
#             return 0.0    
# #        

    else:
            return -np.inf
        
def log_probability(theta, z, mb, model):
    '''returns log probability for a set of values'''
    lp = log_prior(theta, model)
    if not np.isfinite(lp):
        return -np.inf
    if np.isnan(log_likelihood(theta, z, mb, model)):
                print('likelihood is nan!')
    return lp + log_likelihood(theta, z, mb, model)


def run_chains(model, N_step, zarr, m_b):
    '''runs MCMC for given model with N_steps, takes z array, array of mb'''
    
    if model == 'FlatLambdaCDM':
        w_true = -1
        omega_m_true = 0.284
        omega_L_true = .716
        omega_k_true = 0
        initial =np.array([ -19,  50, 0.9])
        pos = initial + 1e-4 * np.random.randn(32, 3) #starting position for chain
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zarr, m_b),  kwargs = {'model':model}
        )
        sampler.run_mcmc(pos, N_step, progress=True);
        samples = (sampler.get_chain(discard= round(N_step/3), flat=True))
        labels = ['M',  'H0', '$\Omega_m$']
        fig = corner.corner(
            flat_samples,  labels=labels, smooth = 1);


        np.savetxt('flat_chains.txt', flat_samples)
        print([np.median(flat_samples[:, 0]),np.median(flat_samples[:, 1]),np.median(flat_samples[:, 2])])
    
   
    elif model == 'wCDM':
        #initial values to start the chain from, ideally use MLE first
        initial =  -19, 60, 0.7, -.9 # M, H0, Om0, w0        
        pos = initial + 1e-4 * np.random.randn(32, 4) #perturb randomly around the initial values
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zCMB, m_b_corr, model)
        )
        sampler.run_mcmc(pos, N_step, progress=True); #run with 50000 steps
        samples = sampler.get_chain(discard= round(N_step/3), thin=15, flat=True) #remove burn in, thin, flatten
        np.savetxt('wCDM_chains_P+.txt', wCDM_samples)
        labels = [  'M', 'H0', '$\Omega_m$',  'w']
        fig = corner.corner(
        wCDM_samples, labels=labels,  smooth = 1);
        
    elif model == 'w0waCDM':
        initial = -19, 70, 0.7, -.9, 0.1 # M, H0, Om0,  w0, wa = theta        
        # pos = soln.x + 1e-4 * np.random.randn(32, )
        pos = initial + 1e-4 * np.random.randn(32, 5)
        nwalkers, ndim = pos.shape

        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zhel, Mb, dMb, model)
        )
        sampler.run_mcmc(pos, 10000, progress=True);
#         wCDM_samples = sample(pos, args, model)

        w0waCDM_samples = sampler.get_chain(discard=500, thin=15, flat=True)
        np.savetxt('w0waCDM_chains.txt', w0waCDM_samples)
        labels = ['M',  'H0', '$\Omega_m$', '$\Omega_0$', 'w0', 'wa']
        fig = corner.corner(
        wCDM_samples, labels=labels,  smooth = 1);
        
    elif model == 'w0wzCDM':   
        initial = -19, 70, 0.7, 0.1, -.9, 0.2 # M, H0, Om0, Ode0, w0, wa = theta        
        # pos = soln.x + 1e-4 * np.random.randn(32, )
        pos = initial + 1e-4 * np.random.randn(32, 6)
        nwalkers, ndim = pos.shape

        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zarr, mb, model)
        )
        sampler.run_mcmc(pos, 10000, progress=True);
#         wCDM_samples = sample(pos, args, model)

        w0wzCDM_samples = sampler.get_chain(discard=500, thin=15, flat=True)
        np.savetxt('w0wzCDM_chains.txt', w0wzCDM_samples)
        labels = ['M',  'H0', '$\Omega_m$', '$\Omega_0$', 'w0', 'wz']
        fig = corner.corner(
        w0wzCDM_samples, labels=labels,  smooth = 1);
        

        
    elif model == 'FLRW_Cosmography':
        #FLRW standard cosmography distance-luminoscity relation
        initial =  70, -19, -.55, .9
#         H0, M, q0, j
        pos = initial + (1e-4 * np.random.randn(32,4))
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zhel, Mb, dMb),  kwargs = {'model':model})
        sampler.run_mcmc(pos, 100000, progress=True);
        flrw_cosmography_samples = (sampler.get_chain(discard=1000, flat=True))
        np.savetxt('FLRW_cosmography_chains.txt', flrw_cosmography_samples)
        labels = ['H0', 'M', 'q0', 'j' ]
        fig = corner.corner(
            flrw_cosmography_samples,  labels=labels, smooth = 1)
        
        
    elif model == 'Cosmography':
        #full cosmography distance-luminoscity relation
        #just Q
#         initial =  70, -19, -.5,0.05, .9, .01
#         pos = initial + (1e-4 * np.random.randn(32,6))
        #   H0_mono, M, q0, j, L1, L2, S= theta #j=Omega_k0 - j0
#         initial =  70, -19, -.5 ,0.5,.1, .1, .01
#         pos = initial + (1e-4 * np.random.randn(32,7))
#quadrupole and dipole

        initial =  70, -19, -.5,0.05, .9, .01, 1., 1., 0.1
#         H0, M, qm, qd, j, S, L1, L2, S2
        pos = initial + (1e-4 * np.random.randn(32,9))
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zhel, Mb, dMb),  kwargs = {'model':model}
        )
        sampler.run_mcmc(pos, 100000, progress=True);
        cosmography_samples = (sampler.get_chain(discard=500, flat=True))
#         np.savetxt('exp_cosmography_chains.txt', cosmography_samples)
#         np.savetxt('cosmography_lin_decay.txt', cosmography_samples)
#          np.savetxt('cosmography_quad.txt', cosmography_samples)
        np.savetxt('H_quad_q_dip.txt', cosmography_samples)
#         labels = ['H0', 'M', 'qm','qd', 'j', 'S' ]
#         labels = ['H0', 'M', 'qm','qm1', 'j', 'S', 'L1', 'L2', 'S2'] 
        labels =  ['H0mono', 'M', 'q0', 'j', 'L1', 'L2', 'S2'] 
        fig = corner.corner(
            cosmography_samples,  labels=labels, smooth = 1)
        
    else:
        raise ValueError
    return(samples)

##################################
#edit here!
##############################


run_chains('wCDM', 100000, zCMB, m_b_corr)
