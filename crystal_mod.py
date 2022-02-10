"""
Author: Wes Johnson
Date:   June 30th 2021
Desc:   This file contains functions developed for veiwing and manipulating 
	the ions crystal in a penning trap."""


#imports: 
import numpy as np
import mode_analysis_code_original
import coldatoms
import ion_trapping
import matplotlib.pyplot as plt
import sys, os
from scipy import constants
from scipy import optimize
import scipy.linalg as LA
import math
import algopy
import os
from numpy import dot 
from numpy import matmul
from scipy.linalg.blas import sgemv
from math import sin,cos
from time import time
import multiprocessing
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter

def checkDir(directory):
    '''Check if directory exists, if not, create it'''
    #make directories if not created: 
    #taken from: 
    #https://djangocentral.com/check-if-a-directory-exists-if-not-create-it/

    # You should change 'test' to your preferred folder.
    MYDIR = (directory)
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)

    else:
        print(MYDIR, "folder already exists.")

def save_ensemble(ddir,num_ions=7,mass_amu=9.012182,
        vtrap=(0.0, -1750.0, -1970.0),v_wall=1.,frot=180000,method='bfgs'):
    """This function generates the initial state of the ion crystal and
        saves it to a file, returning the filename."""
    # Now compute the steady state configuration of the ions for certain trap
    # voltages and rotating wall potential
    frot_kHz = 1e-3*frot 
    mode_analysis=mode_analysis_code_original.ModeAnalysis(
            N      = num_ions,
            ionmass= mass_amu,
            Vtrap  = vtrap,
            Vwall  = v_wall,
            frot   = frot_kHz,
            method = method)
    mode_analysis.run()
    # A few additional parameters
    m_beryllium = mode_analysis.m_Be
    q_beryllium = mode_analysis.q
    # Now create the ensemble
    ensemble = ion_trapping.create_ensemble(mode_analysis.uE,
                                        2.0 * np.pi * frot,
                                        m_beryllium,
                                        q_beryllium)
    # And save it to disk.
    runname = '%dIons' %num_ions
    savename= 'zero_energy_state_'+runname+method+'.txt'
    f = open(ddir+savename, 'w')
    f.write(coldatoms.ensemble_to_json(ensemble))
    f.close()
    return savename

def load_crystal_ensemble(ddir,flnm):
    initial_state = coldatoms.json_to_ensemble(open(ddir+flnm).read())
    return initial_state

def posxxx(posNx3):
    """reshapes array to list all x values, then all y values, then all z values 
        sequentially"""
    return posNx3.copy().T.ravel()

def posxyz(posNx3):
    """reshapes array to list each particles coordinate vector sequentially"""
    return posNx3.copy().reshape((3*N,))

def posNx3(self,posxxx):
    """takes array of x's,y's then z's listed sequentially and turns it into array 
        of particle vector positions or velocities."""
    N = int(np.size(posxxx)/3)
    return np.tile((posxxx[0:N],posxxx[N:2*N],posxxx[2*N:]),reps=1).T

def d(i):
    return i.reshape((i.size, 1)) - i 

def Hii(rsep,rsep5,disq):
    Haiai = np.mat(np.diag(  - np.sum( (rsep ** 2 - 3 * disq) * rsep5, axis=0)))
    Haibi = np.mat(                    (rsep ** 2 - 3 * disq) * rsep5)
    return  Haiai + Haibi

def Hij(rsep5,di,dj):
    Haibi = np.mat(np.diag(  np.sum(3 * di * dj * rsep5, axis=0)))
    Haibj = np.mat(             -   3 * di * dj * rsep5)
    return Haibi + Haibj
   
def kin_rot_eng(velocity):
    shape   = np.shape(velocity)
    indices = len(shape)
    its,ions,dims=shape
    kin_eng_ion = np.empty((its,ions))
    kin_eng_ion = np.sum(np.square(velocity),axis=2)
    return kin_eng_ion

def pot_rot_eng_trap(position,beta,delta):
    """positions must be normalized by ma.l0 from
        mode analysis code"""
    shape   = np.shape(position)
    indices = len(shape)
    its,ions,dims=shape
    x = position[:,:,0]
    y = position[:,:,1]
    z = position[:,:,2]
    pot_eng_ion = np.empty((its,ions))
    pot_eng_ion  = np.square(z)
    pot_eng_ion += beta *(np.square(x)+np.square(y))
    pot_eng_ion += delta*(np.square(x)-np.square(y))
    return pot_eng_ion 
    
def pot_rot_eng_coul(position):
    shape   = np.shape(position)
    indices = len(shape)
    its,ions,dims=shape
    x = position[:,:,0]
    y = position[:,:,1]
    z = position[:,:,2]
    def ds(i):
        return i.reshape((its,ions,1)) - i.reshape((its,1,ions))
    dx,dy,dz, = ds(x),ds(y),ds(z)
    rsep = np.sqrt(dx**2 + dy**2 + dz**2)
    with np.errstate(divide='ignore'):                                                
        rinv = np.where(rsep != 0., rsep ** (-1), 0) 
    pot_eng_ion = np.empty((its,ions))
    pot_eng_ion = np.sum(rinv,axis=1)
    return pot_eng_ion 

def pot_rot_eng(position,beta,delta):
    p_trap = pot_rot_eng_trap(position,beta,delta)
    p_coul = pot_rot_eng_coul(position)
    return p_trap + p_coul 

def coulomb_hessian(pos_array):                                                 
    """Calculate Hessian of potential"""                                              
    Nion = int(np.size(pos_array)/3)
    #print('nion: %d' %Nion)
    l0=1e-5
    C = constants.e**2 / (4*np.pi*constants.epsilon_0)/l0**3
    x = pos_array[0:Nion]/l0                                                        
    y = pos_array[Nion:2*Nion]/l0
    z = pos_array[2*Nion:]/l0

    dx = d(x)                                                 
    dy = d(y) 
    dz = d(z)#for planar crystal this is always zero
    
    #print(dx,'\n')
    #print(dy,'\n')    
    #print(dz,'\n')        
    
    dxsq = dx ** 2                                                                    
    dysq = dy ** 2                           
    dzsq = dz ** 2 #for planar crystal this is zero
    #print(dxsq,'\n')
    #print(dysq,'\n')    
    #print(dzsq,'\n')        
    rsep = np.sqrt(dxsq + dysq + dzsq) 
    #print(rsep,'\n') 
    with np.errstate(divide='ignore'):                                                
        rsep5 = np.where(rsep != 0., rsep ** (-5), 0) 
    #Double derivatives   
    Hxx = Hii(rsep,rsep5,dxsq)                                      
    Hyy = Hii(rsep,rsep5,dysq)                                      
    Hzz = Hii(rsep,rsep5,dzsq) 
    #print(Hxx,'\n')
    #print(Hyy,'\n')
    #print(Hzz,'\n')
    

    #Mixed derivatives                                                               
    Hxy = Hij(rsep5,dx,dy)
    Hyx = Hxy
    Hxz = Hij(rsep5,dx,dz) 
    Hzx = Hxz
    Hyz = Hij(rsep5,dy,dz)
    Hzy = Hyz
    #print(Hxy,'\n')
    #print(Hxz,'\n')
    #print(Hyz,'\n')
    

    H = np.bmat([[Hxx, Hxy, Hxz], [Hyx, Hyy, Hyz], [Hzx, Hzy, Hzz]])
    H = np.asarray(H)*C
    return H

class LinearizedCoulombForce(object):
    """This class produces a force from a force matrix based on derivatives of the 
        coulomb potential"""

    def __init__(self,wrot,matrix,ensemble):
        """For the planar crystal the equilibrium position in z is just zero so we only
            need the hessian matrix."""
        self.N         = ensemble.num_ptcls
        self.pos0      = np.array(ensemble.x.copy())
        self.pos       = np.empty_like(ensemble.x)
        self.rot_vec   = np.empty_like(ensemble.x)
        self.pos_temp  = np.empty_like(ensemble.x)
        self.rot_pos   = np.empty_like(ensemble.x)

        self.time      = 0.0
        self.iteration = 0
        self.omega_rot = wrot
        self.rot_mat   = np.empty((2,2))

        self.force_mat = matrix
        self.force_vec = np.zeros_like(ensemble.x)
        #use coldatoms library to get the force at equilibrium
        coldatoms.CoulombForce().force(1,ensemble,self.force_vec)

    def posxxx(self,posNx3):
        """reshapes array to list all x values, then all y values, then all z values 
            sequentially"""
        return posNx3.T.ravel()

    def posNx3(self,posxxx):
        """takes array of x's,y's then z's listed sequentially and turns it into array 
            of particle vector positions or velocities."""
        return posxxx.reshape(self.N,3,order='F')

    def get_rot_mat(self,theta):
        c,s = cos(theta),sin(theta)
        self.rot_mat = np.array([[c,-s],[s,c]])

    def rotateNx3(self,pos,rot_mat):
        self.pos_temp[:,0:2] = np.matmul(pos[:,0:2],rot_mat)
        self.pos_temp[:,2]   = pos[:,2]
        return self.pos_temp

    def get_rot_pos(self):
        self.rot_pos[:,0:2] = np.matmul(self.pos[:,0:2],self.rot_mat)
        self.rot_pos[:,2]   = self.pos[:,2]

    def force(self, dt, ensemble, f):

        self.pos = ensemble.x

        self.iteration += 1
        self.time = float(self.iteration)*dt

        theta          = self.time*self.omega_rot
        self.get_rot_mat(theta)
        self.get_rot_pos()
        #self.rot_pos   = self.rotateNx3(self.pos,self.rot_mat)
        self.rot_vec   = self.posxxx(self.rot_pos - self.pos0)
        f0Plusf1       = self.rotateNx3(self.force_vec - 
                            self.posNx3(matmul(self.force_mat,self.rot_vec)),self.rot_mat.T )

        f += (dt)*f0Plusf1



class LinearizedCoulombForce_copy(object):
    """This class produces a force from a force matrix based on derivatives of the 
        coulomb potential"""

    def __init__(self,wrot,matrix,ensemble):
        """For the planar crystal the equilibrium position in z is just zero so we only
            need the hessian matrix."""

        self.N         = ensemble.num_ptcls
        self.pos       = np.array(ensemble.x.copy())

        self.time      = 0.0
        self.iteration = 0
        self.omega_rot = wrot

        self.force_mat = matrix
        self.force_vec = np.zeros_like(ensemble.x)
        #use coldatoms library to get the force
        coldatoms.CoulombForce().force(1,ensemble,self.force_vec)

    def posxxx(self,posNx3):
        """reshapes array to list all x values, then all y values, then all z values 
            sequentially"""
        return posNx3.copy().T.ravel()

    def posxyz(self,posNx3):
        """reshapes array to list each particles coordinate vector sequentially"""
        return posNx3.copy().reshape((3*N,))

    def posNx3(self,posxxx):
        """takes array of x's,y's then z's listed sequentially and turns it into array 
            of particle vector positions or velocities."""
        N = int(np.size(posxxx)/3)
        return np.tile((posxxx[0:N],posxxx[N:2*N],posxxx[2*N:]),reps=1).T

    def get_rot_mat(self,theta):
        return np.array([[cos(theta),-sin(theta)],[sin(theta),cos(theta)]])

    def rotateNx3(self,pos,rot_mat):
        rot_pos = np.empty((self.N,3))
        rot_pos[:,0:2] = np.matmul(pos[:,0:2],rot_mat)
        rot_pos[:,2]   = pos[:,2]
        return rot_pos

    def force(self, dt, ensemble, f):

        pos = ensemble.x

        self.iteration += 1
        self.time = float(self.iteration)*dt

        theta          = self.time*self.omega_rot
        rot_mat        = self.get_rot_mat(theta)

        rot_pos        = self.rotateNx3(pos,rot_mat)
        rot_vec        = self.posxxx(rot_pos - self.pos)
        f0Plusf1       = self.rotateNx3(self.force_vec - 
                            self.posNx3(matmul(self.force_mat,rot_vec)),rot_mat.T )

        f += (dt)*f0Plusf1



class TrapPotential(object):

    def __init__(self, kz, delta, omega, phi_0):
        self.kz = kz
        self.kx = -(0.5 + delta) * kz
        self.ky = -(0.5 - delta) * kz
        self.phi_0 = phi_0
        self.phi = phi_0
        self.omega = omega

    def reset_phase(self):
        self.phi = self.phi_0
            
    def force(self, dt, ensemble, f):
        self.phi += self.omega * 0.5 * dt
        
        q = ensemble.ensemble_properties['charge']
        if q is None:
            q = ensemble.particle_properties['charge']
            if q is None:
                raise RuntimeError('Must provide ensemble or per particle charge')

        cphi = np.cos(self.phi)
        sphi = np.sin(self.phi)
        kx = self.kx
        ky = self.ky
        
        x = ensemble.x[:, 0]
        y = ensemble.x[:, 1]
        z = ensemble.x[:, 2]
        
        f[:, 0] += dt * q * (
            (-kx * cphi * cphi - ky * sphi * sphi) * x + cphi * sphi * (ky - kx) * y)
        f[:, 1] += dt * q * (
            cphi * sphi * (ky - kx) * x + (-kx * sphi * sphi - ky * cphi * cphi) * y)
        f[:, 2] += -dt * q *self.kz * z

        self.phi += self.omega * 0.5 * dt



def evolve_ensemble(dt, t_max, ensemble, Bz, forces):
    num_steps = int(t_max / dt)-1
    coldatoms.bend_kick(dt, Bz, ensemble, forces, num_steps=num_steps)
    coldatoms.bend_kick(t_max - dt * num_steps, Bz, ensemble, forces)

def integrateTrajectories(iterations,dt, t_max, ensemble, B,forces):
    start = time()   
    N = ensemble.num_ptcls
    times      = np.empty((iterations+1))
    times[0]   = 0.0
    posall     = np.empty((iterations+1,N,3))
    posall[0,:,:]=ensemble.x[:,:] 
    for it in range(1,iterations+1):
        evolve_ensemble(dt, t_max, ensemble, B, forces)
        posall[it,:,:]       = ensemble.x[:,:]
        times[it]            = it*t_max 
    end = time()
    duration = end - start 
    print("Time: %1.4f"%duration)
    return posall


def integrateStates(steps, dt, t_max, ensemble, B,forces,save_step=1):
    iterations = (steps +1)* save_step 
    N          = ensemble.num_ptcls
    step       = 0
    state_type = np.dtype([('position', np.float64,  (steps+1,N,3)), 
                            ('velocity', np.float64, (steps+1,N,3)),
                            ('time',np.float64,      (steps+1,)
                            )])
    start = time()   
    times      = np.empty((steps+1))
    times[0]   = 0.0
    posall     = np.empty((steps+1,N,3))
    velall     = np.empty((steps+1,N,3))
    posall[0,:,:]=ensemble.x[:,:] 
    velall[0,:,:]=ensemble.v[:,:] 
    for it in range(1,iterations):
        evolve_ensemble(dt, t_max, ensemble, B, forces)
        if ((it)%save_step==0):
            step += 1 
            posall[step,:,:]       = ensemble.x[:,:]
            velall[step,:,:]       = ensemble.v[:,:]
            times[step]            = it*t_max 
            #print('step = %d, it = %d'%(step,it))
    end = time()
    duration = end - start 
    print("Time: %1.4f"%duration)
    return np.array([(posall,velall,times)], dtype=state_type)

def integrateStates_copy(iterations, dt, t_max, ensemble, B,forces):
    N = ensemble.num_ptcls
    state_type = np.dtype([('position', np.float64,  (iterations+1,N,3)), 
                            ('velocity', np.float64, (iterations+1,N,3)),
                            ('time',np.float64,      (iterations+1,)
                            )])
    start = time()   
    times      = np.empty((iterations+1))
    times[0]   = 0.0
    posall     = np.empty((iterations+1,N,3))
    velall     = np.empty((iterations+1,N,3))
    posall[0,:,:]=ensemble.x[:,:] 
    velall[0,:,:]=ensemble.v[:,:] 
    for it in range(1,iterations+1):
        evolve_ensemble(dt, t_max, ensemble, B, forces)
        posall[it,:,:]       = ensemble.x[:,:]
        velall[it,:,:]       = ensemble.v[:,:]
        times[it]            = it*t_max 
    end = time()
    duration = end - start 
    print("Time: %1.4f"%duration)
    return np.array([(posall,velall,times)], dtype=state_type)


def psd(coord_traj,time_step,units=1e6,log=True):
    """computes the power spectrum of the coordianate trajectroy and converts the units of the frequecies from hertz
        default is MHz. Returns the log of the power spectrum by default, and the frequencies."""
    ps =  np.sum(np.abs(np.fft.fft(coord_traj,axis=0) / np.size(coord_traj[0,:]))**2, axis=1)
    freqs = np.fft.fftfreq(coord_traj[:,0].size, time_step)
    idx = np.argsort(freqs)
    freqs_norm=freqs[idx]/units
    if log == True:
        ps = np.log(ps[idx])
    return ps,freqs_norm




def rotating_frame(posall,omega=2*np.pi*180e3,t_max=1e-8):
    posall_rot=np.empty_like(posall)
    thetas = -np.arange(np.size(posall[:,0,0]))*t_max*omega
    N = np.size(posall[0,:,0])
    cs1    = np.cos(thetas)
    sn1    = np.sin(thetas)
    cs     = np.tile(cs1,(N,1)).T
    sn     =  np.tile(sn1,(N,1)).T
    posall_rot[:,:,0] = cs*posall[:,:,0] - sn*posall[:,:,1]
    posall_rot[:,:,1] = sn*posall[:,:,0] + cs*posall[:,:,1]
    return posall_rot


def randomize_phase(evects):
    N = int(np.size(evects[0,:])/4)
    phase = 2*np.pi*np.random.random(N)-np.pi
    complexamps = np.exp(complex(0,1)*phase)

    randomphases = np.empty_like(evects)

    for vector in np.arange(0,N,1):
        randomphases[:,vector] = evects[:,vector]*complexamps[vector]
    return randomphases


def initialize_planar_modes(ensemble,evects,amps):
    new_ensemble = ensemble.copy()
    N = new_ensemble.num_ptcls
    x = new_ensemble.x[:,0]
    y = new_ensemble.x[:,1]
    vx= new_ensemble.v[:,0]
    vy= new_ensemble.v[:,1]
    for i in range(2*N):#there are 2*N planar modes
        ev = evects[:,2*i]#the 2*i+1 vector is c.c.
        new_ensemble.x[:,0] = x + amps[i]*np.real(ev[0  :  N].T)
        new_ensemble.x[:,1] = y + amps[i]*np.real(ev[  N:2*N].T)
        new_ensemble.v[:,0] = vx+ amps[i]*np.real(ev[2*N:3*N].T)
        new_ensemble.v[:,1] = vy+ amps[i]*np.real(ev[3*N:4*N].T)
    return new_ensemble
#ccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cc
#laser cooling functions 
class UniformBeam(object):
    """A laser beam with a uniform intensity profile."""

    def __init__(self, S0):
        """Construct a Gaussian laser beam from position, direction, and width.

        S0 -- Peak intensity (in units of the saturation intensity).
        k -- Propagation direction of the beam (need not be normalized).
        """
        self.S0 = S0

    def intensities(self, x):
        return self.S0 * np.ones_like(x[:, 0])

class DopplerDetuning(object):
    def __init__(self, Delta0, k):
        self.Delta0 = Delta0
        self.k = np.copy(k)

    def detunings(self, x, v):
        return self.Delta0 - np.inner(self.k, v)

class GaussianBeam(object):
    """A laser beam with a Gaussian intensity profile."""

    def __init__(self, S0, x0, k, sigma):
        """Construct a Gaussian laser beam from position, direction, and width.

        S0 -- Peak intensity (in units of the saturation intensity).
        x0 -- A location on the center of the beam.
        k -- Propagation direction of the beam (need not be normalized).
        sigma -- 1/e width of the beam."""
        self.S0 = S0
        self.x0 = np.copy(x0)
        self.k_hat = k / np.linalg.norm(k)
        self.sigma = sigma

    def intensities(self, x):
        xp = x - self.x0
        xperp = xp - np.outer(xp.dot(self.k_hat[:, np.newaxis]), self.k_hat)
        return self.S0 * np.exp(-np.linalg.norm(xperp, axis=1)**2 / self.sigma**2)

#ccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cc
# Eigen-mode diagnostics/initialization 
# initialization: 

def initialize_planar_eigenmode(ensemble,ma,evs,modes,amp,wrot):
    """The function initializes the planar eigenmodes of the given
        ensemble.
        Variables:
        ensemble :: cold atoms object 
        ma       :: instance of mode analysis code
        evs      :: the eigenvectors (with or without
                    random phase)
        modes    :: list of modes to be excited
        amp      :: amplitude of the modes
    """
    eigen_state = ensemble.copy()    
    N           = eigen_state.num_ptcls
    for mode in modes: 
        
        dx = np.real(evs[:,2*mode][0*N:1*N].T)*amp*ma.l0
        dy = np.real(evs[:,2*mode][1*N:2*N].T)*amp*ma.l0
        dvx= np.real(evs[:,2*mode][2*N:3*N].T)*amp*ma.v0
        dvy= np.real(evs[:,2*mode][3*N:4*N].T)*amp*ma.v0

        
        x  = eigen_state.x[:,0]
        y  = eigen_state.x[:,1]
        vx = eigen_state.v[:,0]
        vy = eigen_state.v[:,1]


        eigen_state.x[:,0] = x + dx 
        eigen_state.x[:,1] = y + dy
        eigen_state.v[:,0]= vx +  dvx - wrot*dy 
        eigen_state.v[:,1]= vy +  dvy + wrot*dx

    return eigen_state

def initialize_axial_eigenmode(ensemble,ma,evs,modes,amp):
    """The function initializes the axial eigenmodes of the given
        ensemble.
        Variables:
        ensemble :: cold atoms object 
        ma       :: instance of mode analysis code
        evs      :: the eigenvectors (with or without
                    random phase)
        modes    :: list of modes to be excited
        amp      :: amplitude of the modes
    """
    eigen_state = ensemble.copy()
    N           = eigen_state.num_ptcls
    for mode in modes: 
        
        dz = np.real(evs[:,2*mode][0*N:1*N].T)*amp*ma.l0
        dvz= np.real(evs[:,2*mode][1*N:2*N].T)*amp*ma.v0
        
        z  = eigen_state.x[:,2]
        vz = eigen_state.v[:,2]

        eigen_state.x[:,2] = z + dz  
        eigen_state.v[:,2] =vz + dvz
        
    return eigen_state


def randomize_phases(Evects):
    """Randomizes the complex phases of the given array of 
        eigenvectors, returns the phases and the array of 
        randomized phase eigenvectors. 
        Variables: 
        Evects      :: array of eigen vectors
    """
    Pt     = Evects.copy()
    N      = np.size(Pt[0,:])
    thetas = 2*np.pi*np.random.random((N,))
    phases = np.exp(1j*thetas)
    for i in range(N):
        Pt[:,i] *= phases[i]
    return Pt,phases

def temp2amp(temperature,ma):
    """Given mode analysis and a temperature, the function returns the 
        corresponding amplitude. 
        Variables: 
        temperature     :: temperture of mode in milikelvin 
        ma              :: mode analysis instance 
    """
    return np.sqrt(constants.k * 2*temperature / ma.E0)

def amp2temp(amp,ma):
    """Given mode analysis and an amplitude, the function returns the 
        corresponding temperture Kelvin. 
        Variables: 
        amp             :: mode analysis units of amplitude  
        ma              :: mode analysis instance 
    """
    return 1/2*np.square(amp)*ma.E0/constants.k
    
def get_axial_state(state,mode_analysis):
    ma          = mode_analysis
    stsq        = np.squeeze(state)
    posz        = stsq['position'][:,:,2]
    velz        = stsq['velocity'][:,:,2]
    axial_state = np.concatenate((
                                posz/ma.l0,
                                velz/ma.v0),
                                axis = -1)
    return axial_state

def get_planar_state(state,equ_pos,mode_analysis,omega,dt):
    ma        = mode_analysis
    stsq      = np.squeeze(state)
    posall    = stsq['position']  
    velall    = stsq['velocity'] 

    posall_rot= np.empty_like(posall)
    velall_rot= np.empty_like(velall)
    thetas    = -np.arange(np.size(posall[:,0,0]))*dt*omega
    N         = np.size(posall[0,:,0])
    cs1       = np.cos(thetas)
    sn1       = np.sin(thetas)
    cs        = np.tile(cs1,(N,1)).T
    sn        = np.tile(sn1,(N,1)).T
    
    vel_col_rot = np.empty_like(velall)
    vel_col_rot[:,:,0] = velall[:,:,0] + omega*posall[:,:,1]
    vel_col_rot[:,:,1] = velall[:,:,1] - omega*posall[:,:,0]
    
    posall_rot[:,:,0] = cs*posall[:,:,0] - sn*posall[:,:,1]
    posall_rot[:,:,1] = sn*posall[:,:,0] + cs*posall[:,:,1]
    
    velall_rot[:,:,0] = cs*vel_col_rot[:,:,0] - sn*vel_col_rot[:,:,1]
    velall_rot[:,:,1] = sn*vel_col_rot[:,:,0] + cs*vel_col_rot[:,:,1]
    
    planar_state = np.concatenate(
                                (
                                (posall_rot[:,:,0] - equ_pos[:,0])/ma.l0,
                                (posall_rot[:,:,1] - equ_pos[:,1])/ma.l0,
                                (velall_rot[:,:,0])/ma.v0,
                                (velall_rot[:,:,1])/ma.v0                                    
                                ),
                                axis = -1)
    return planar_state

def get_eig_amp(eng_mat,eig_vec,state):
    eVs     = np.conj(eig_vec).T
    Z       = state.T
    E_Z     = np.matmul(eng_mat,Z)
    eVs_E_Z = np.matmul(eVs,E_Z)
    return 2*np.abs(eVs_E_Z).T

def get_eig_amp_conjplusvec(eng_mat,eig_vec,state):
    m       = np.arange(0,np.size(state[0,:]),2)
    eVs     = np.conj(eig_vec).T
    Z       = state.T
    E_Z     = np.matmul(eng_mat,Z)
    eVs_E_Z = np.matmul(eVs,E_Z)
    vamps   = (eVs_E_Z)
    mamps   = vamps.copy()
    mamps[m]  = vamps[m] + vamps[m+1]
    mamps[m+1]= vamps[m] + vamps[m+1]
    return np.real(mamps).T


def get_axial_eng_mat(ma):
    return LA.block_diag(ma.axial_hessian,ma.axial_Mmat)

def get_planar_eng_mat(ma):
    return LA.block_diag(ma.planar_hessian,ma.planar_Mmat)


def sort_eig_amps(amps):
    sq_amps     = np.squeeze(amps)
    num_amps    = np.shape(sq_amps)[1]
    unique_amps = (sq_amps)[:,np.arange(0,num_amps,2)]
    amp_means   = np.mean(unique_amps,axis = 0)
    sort_amps   = np.ravel(np.argsort(-amp_means))
    return sort_amps

def sort_amps_std(amps): 
    sq_amps     = np.squeeze(amps)
    num_amps    = np.shape(sq_amps)[1]
    unique_amps = (sq_amps)[:,np.arange(0,num_amps,2)]
    amp_stds    = np.std(unique_amps,axis = 0)
    sort_amps   = np.ravel(np.argsort(-amp_stds))
    return sort_amps


#ccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cc
#plotting functions 
def plot_crystal_xy(pos):
    x = pos[:,0]
    y = pos[:,1]
    fig = plt.figure()
    plt.scatter(x=x/1e-6,y=y/1e-6)
    plt.title("Planar Veiw of Ion Crystal")
    plt.xlabel(r"x [$\mu$m]")
    plt.ylabel(r"y [$\mu$m]")
    lim = 1.1*np.max(np.array(pos[:,0:2]))/1e-6
    plt.axis([-lim, lim, -lim, lim])
    plt.axis('equal')

def plot_crystal_xy_ax(pos,ax,lim=None):
    x = pos[:,0]
    y = pos[:,1]
    ax.scatter(x=x/1e-6,y=y/1e-6
            ,label='Positions')
    ax.set_title("Planar Veiw of Ion Crystal")
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel(r"y [$\mu$m]")
    if lim==None:
        lim = 1.1*np.max(np.array(pos[:,0:2]))/1e-6
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    ax.set_aspect('equal', 'box')

def plot_state_lab(state,ax,fig):
    x = state.x[:,0]*1e6
    y = state.x[:,1]*1e6
    N = state.num_ptcls
    lim = np.max([x,y])*1.2
    for i in range(N):
        dvx = state.v[i,0]/N*3
        dvy = state.v[i,1]/N*3
        ax.arrow(x[i],y[i],dvx,dvy,width=1,fc='black')
    ax.set_aspect(1)
    ax.scatter(x,y,label='$\mathbf{X}$')
    ax.set_xlabel('x [$\mu m$]')
    ax.set_ylabel('y [$\mu m$]')
    ax.set_title('State in Lab Frame')
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.scatter(10**4,10**4,marker=r'$\rightarrow$',color='k'
              ,label=r'$\mathbf{V}$',s=250)
    ax.legend(loc='upper left')

def get_mode_ticks(labels,freqs,min_diff=0.05):
    assert len(labels)==len(freqs),"labels and freqs must be same length"
    diffs = np.diff(freqs)/(freqs[-1]-freqs[0])
    gdexs = np.append(np.where(diffs>min_diff),np.array([len(freqs)-1]))
    bdexs = np.array(np.where(diffs<=min_diff))
    tdexs = np.isin(gdexs,bdexs+1)
    fdexs = gdexs[np.where(np.logical_not(tdexs))]
    stops = gdexs[np.where(tdexs)]
    bgins = bdexs[np.where(np.isin(bdexs,gdexs+1))]
    assert len(bgins)==len(stops),'something went wrong'
    dexs  = np.append(gdexs,stops)
    freqs[stops] = (freqs[stops] + freqs[bgins])/2 
    labels[stops]  = [labels[bgins[i]]+'-'+labels[stops[i]] for i in range(len(bgins))]
    gfreqs = freqs[dexs] ; glabels = labels[dexs]
    return gfreqs,glabels 

def plot_planar_eigenmode(mode_analysis,mode,ax,fig,theta=0):
    ma            = mode_analysis
    N             = ma.Nion
    posxxx        = np.append(ma.uE,np.repeat([0.0],N))
    posNx3        = posxxx.reshape((N,3),order = 'F')
    equ_pos       = posNx3
    x             = posNx3[:,0]*1e6
    y             = posNx3[:,1]*1e6

    evs = ma.planarEvects
    om  = ma.planarEvalsE 
    ev = -evs[:,mode*2]*np.exp(complex(0,theta))

    ax.scatter(x=x,y=y,color='royalblue',zorder = 3)
    ax.set_aspect('equal', 'box')
    ax.set_title("n = %d, "%(mode+1)+r"$\omega_n$"+"=%1.2e[Hz]" %(om[mode]))
    lim = np.max(np.abs([x,y]))*1.25
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)   
    dx  = np.real(ev)[0*N:1*N]
    dy  = np.real(ev)[1*N:2*N]
    dvx = np.real(ev)[2*N:3*N]
    dvy = np.real(ev)[3*N:4*N]
    norm = np.linalg.norm(   np.concatenate(  (dx,dy,dvx,dvy)   )  )/(2*N/2.5)
    if norm>0:
        dx /= norm 
        dy /= norm 
        dvx/= norm 
        dvy/= norm 
    for ion in range(N):
        ax.arrow(x[ion],y[ion],dx[ion,0],dy[ion,0],width = 1,head_length=2,fc='black', ec='black')
        ax.arrow(x[ion],y[ion],dvx[ion,0],dvy[ion,0],width = 1,head_length=2,fc='red', ec='red')

    fig.set_size_inches(8, 8)
    leg = fig.legend([r'$\mathbf{\delta x}$'
                        ,r'$\mathbf{\delta v}$'
                        ,r'$\mathbf{X_0}$']
                        ,labelcolor = ['black','red','royalblue'], 
                     loc = 'upper left')
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('red')
    leg.legendHandles[2].set_color('royalblue')
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel(r"y [$\mu$m]")

def plot_planar_eigenmode_pos(initial_state,mode_analysis,mode,ax,fig,theta=0):
    ma = mode_analysis
    x = initial_state[:,0]*1e6
    y = initial_state[:,1]*1e6
    N = np.size(initial_state[:,0])

    evs = ma.planarEvects
    om  = ma.planarEvalsE 
    ev = -evs[:,mode*2]*np.exp(complex(0,theta))

    ax.scatter(x=x,y=y,color='royalblue',zorder = 3)
    ax.set_aspect('equal', 'box')
    ax.set_title("n = %d, "%(mode+1)+r"$\omega_n$"+"=%1.2e[Hz]" %(om[mode]))
    lim = np.max(np.abs([x,y]))*1.25
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)   
    dx  = np.real(ev)[0*N:1*N]
    dy  = np.real(ev)[1*N:2*N]
    dvx = np.real(ev)[2*N:3*N]
    dvy = np.real(ev)[3*N:4*N]
    norm = np.linalg.norm(   np.concatenate(  (dx,dy,dvx,dvy)   )  )/(2*N/2.5)
    if norm>0:
        dx /= norm 
        dy /= norm 
        dvx/= norm 
        dvy/= norm 
    for ion in range(N):
        ax.arrow(x[ion],y[ion],dx[ion,0],dy[ion,0],width = 1,head_length=2,fc='black', ec='black')
        ax.arrow(x[ion],y[ion],dvx[ion,0],dvy[ion,0],width = 1,head_length=2,fc='red', ec='red')

    fig.set_size_inches(8, 8)
    leg = fig.legend([r'$\mathbf{\delta x}$'
                        ,r'$\mathbf{\delta v}$'
                        ,r'$\mathbf{X_0}$']
                        ,labelcolor = ['black','red','royalblue'], 
                     loc = 'upper left')
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('red')
    leg.legendHandles[2].set_color('royalblue')
    ax.set_xlabel(r"x [$\mu$m]")
    ax.set_ylabel(r"y [$\mu$m]")

def plot_psd(freqs_norm,ps_log,ax,fig,
        title=None,label=None,normal_freqs=[],
        color='b'):
    if label==None: 
        label = 'PSD'
    if title ==None:
        title = 'Power Spectrum'
    ax.plot(freqs_norm, ps_log,label=label,color=color)
    if len(normal_freqs)!=0:
        switch = 0
        for norm_freq in normal_freqs:
            if switch == 0: 
                ax.axvline(x=norm_freq,color='r',label='Linear Modes'
                           ,alpha = 0.75,zorder = 1)
                switch += 1
            else: 
                ax.axvline(x=norm_freq,color='r',alpha = 0.75
                           ,zorder = 1)
        nf_min = np.min(normal_freqs);nf_max=np.max(normal_freqs)
        expand = .05*(nf_max - nf_min)
        ax.set_xlim(( nf_min-expand , nf_max+expand ))
    ax.set_ylim((-30, 0))
    ax.set_ylabel(r"$\ln(PSD)$ [A.U.]")
    ax.set_xlabel(r"Frequency [$MHz$]")
    ax.set_title(title)
    

class animate_nlVlin_rotating_frame:
    def __init__(self,posall_nl,posall,wrot,dt,save=False):
        max_slices = 5*10**5
        len_slices = np.size(posall[:,0,0])
        if len_slices<max_slices: 
            sl=len_slices
        else: 
            sl=max_slices
        
        self.posall_rot_nl = rotating_frame(posall_nl[:sl,:,:],
                                            omega=wrot,
                                            t_max=dt)
        self.posall_rot    = rotating_frame(posall[:sl,:,:],
                                            omega=wrot,
                                            t_max=dt)


        self.x = []
        self.y = []

        self.positions1 = self.posall_rot[:,:,0:2]*1e6
        self.positions2 = self.posall_rot_nl[:,:,0:2]*1e6


        self.figure, self.ax = plt.subplots(figsize=(7,7))
        self.lim = np.max(np.abs(posall[0,:,:]*1e6))*1.2
        self.points1, = self.ax.plot(self.x, self.y
                                    ,'ro'
                                    ,label='Linear'
                                    ,markersize=5)
        self.points2, = self.ax.plot(self.x, self.y
                                    ,'bo'
                                    ,label='Nonlinear'
                                    ,markersize=5)
        self.lines1, = self.ax.plot(self.x, self.y,'ro',markersize=0.1)
        self.lines2, = self.ax.plot(self.x, self.y,'bo',markersize=.1)
        self.ax.legend()
        plt.title("Comparison in Rotating Frame")
        self.ax.set_xlabel(r'x [$\mu m$]')
        self.ax.set_ylabel(r'y [$\mu m$]')

    def func_animate(self,frame_number,skip=200,tail=10**5,dots_every=100):
        self.ax.set_ylim((-self.lim,self.lim))
        self.ax.set_xlim((-self.lim,self.lim))
        x1 = self.positions1[frame_number*skip,:,0]
        y1 = self.positions1[frame_number*skip,:,1]
        x2 = self.positions2[frame_number*skip,:,0]
        y2 = self.positions2[frame_number*skip,:,1]
        if frame_number*skip-tail >= 0:
            index = np.arange(frame_number*skip-tail,
                                frame_number*skip,
                                dots_every)
            linex1=self.positions1[index,:,0]
            liney1=self.positions1[index,:,1]
            linex2=self.positions2[index,:,0]
            liney2=self.positions2[index,:,1]
        else:
            index = np.arange(0,frame_number*skip,dots_every)
            linex1=self.positions1[index,:,0]
            liney1=self.positions1[index,:,1]
            linex2=self.positions2[index,:,0]
            liney2=self.positions2[index,:,1]

        thetime = 1e-8*float(frame_number*skip)/1e-6

        self.points1.set_data(x1,y1)
        self.points2.set_data(x2,y2)
        self.lines1.set_data(linex1,liney1)
        self.lines2.set_data(linex2,liney2)

class animate_rotating_frame:
    def __init__(self,posall,wrot,dt,units=1e6):
        max_slices = 10**5
        len_slices = np.size(posall[:,0,0])
        if len_slices<max_slices: 
            sl=len_slices
        else: 
            sl=max_slices
        
        self.posall_rot    = rotating_frame(posall[:sl,:,:],
                                            omega=wrot,
                                            t_max=dt)


        self.x = []
        self.y = []

        self.positions1 = self.posall_rot[:,:,0:2]*units


        self.figure, self.ax = plt.subplots(figsize=(16,16))
        self.lim = np.max(np.abs(posall[0,:,:]*units))*1.2
        self.points1, = self.ax.plot(self.x, self.y
                                    ,'bo'
                                    ,zorder=2
                                    ,label='Particle'
                                    ,markersize=5)
        self.lines1, = self.ax.plot(self.x, self.y
                ,'ro'
                ,zorder=1
                ,label='Trail'
                ,markersize=0.1)
        self.ax.legend()
        plt.title("Rotating Frame Movie")
        self.ax.set_xlabel(r'x [$\mu m$]')
        self.ax.set_ylabel(r'y [$\mu m$]')
        self.ax.set_aspect(1)

    def func_animate(self,frame_number,skip=10,tail=10**3,dots_every=1):
        self.ax.set_ylim((-self.lim,self.lim))
        self.ax.set_xlim((-self.lim,self.lim))
        x1 = self.positions1[frame_number*skip,:,0]
        y1 = self.positions1[frame_number*skip,:,1]
        if frame_number*skip-tail >= 0:
            index = np.arange(frame_number*skip-tail,
                                frame_number*skip,
                                dots_every)
            linex1=self.positions1[index,:,0]
            liney1=self.positions1[index,:,1]
        else:
            index = np.arange(0,frame_number*skip,dots_every)
            linex1=self.positions1[index,:,0]
            liney1=self.positions1[index,:,1]

        thetime = 1e-8*float(frame_number*skip)/1e-6
        #self.ax.annotate(r"Time =%1.0f[$\mu s$]" %thetime, xy=(0.1, 0.04),\
        #            xycoords='axes fraction',size=10,\
        #            bbox=dict(boxstyle='round', fc='w'))
        self.points1.set_data(x1,y1)
        self.lines1.set_data(linex1,liney1)

#ccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cccc5ccc10cc






if __name__ == "__main__":
    pass 
