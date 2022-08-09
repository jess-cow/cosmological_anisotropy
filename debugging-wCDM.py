from astropy.cosmology import wCDM
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
N_step = 100000
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


def log_likelihood(theta, z, m_b_corr,model):
    '''log likelihood function m '''
    if model == 'wCDM' :
        M, H0, Om0, w0 = theta
        Ode0 = 1 - Om0
        flc = wCDM(H0, Om0, Ode0, w0)
        mod = flc.distmod(z).value + M 

    delta_m =m_b_corr - mod
    chi2 = np.matmul(np.transpose(delta_m), np.matmul(invC, delta_m)) #equation8
    L = -0.5*chi2

    
    return(L)

def log_prior(theta, model):
    '''priors on sampled parameters for each model, returns -inf if outside of priors, 0 if inside'''

    if model =='wCDM':
        M, H0, Om0, w0, = theta
        if 0.0 <= Om0 < 1.0  and -20. < M  < 20. and  20. < H0 < 100. and -5 < w0 < 0.:
            return 0.0
        else :
            return -np.inf

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
   
    if model == 'wCDM':
        #initial values to start the chain from, ideally use MLE first
        initial =  -19, 60, 0.7, -.9 # M, H0, Om0, w0        
        pos = initial + 1e-4 * np.random.randn(32, 4) #perturb randomly around the initial values
        nwalkers, ndim = pos.shape
        sampler = emcee.EnsembleSampler(
            nwalkers, ndim, log_probability, args=(zCMB, m_b_corr, model)
        )
        sampler.run_mcmc(pos, N_step, progress=True); #run with 50000 steps
        wCDM_samples = sampler.get_chain(discard= round(N_step/3), thin=15, flat=True) #remove burn in, thin, flatten
        np.savetxt('wCDM_chains_P+.txt', wCDM_samples)
        labels = [  'M', 'H0', '$\Omega_m$',  'w']
        fig = corner.corner(
        wCDM_samples, labels=labels,  smooth = 1);
        
    else:
        raise ValueError
    return(samples)


model = 'wCDM'
initial =  -19, 60, 0.7, -.9 # M, H0, Om0, w0        
pos = initial + 1e-4 * np.random.randn(32, 4) #perturb randomly around the initial values
nwalkers, ndim = pos.shape
sampler = emcee.EnsembleSampler(
nwalkers, ndim, log_probability, args=(zCMB, m_b_corr, model)
)
sampler.run_mcmc(pos, N_step, progress=True); #run with 50000 steps
wCDM_samples = sampler.get_chain(discard= round(N_step/3), thin=15, flat=True) #remove burn in, thin, flatten
np.savetxt('wCDM_chains_P+.txt', wCDM_samples)
labels = [  'M', 'H0', '$\Omega_m$',  'w']
fig = corner.corner(
wCDM_samples, labels=labels);

