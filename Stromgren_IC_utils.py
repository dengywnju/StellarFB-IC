import numpy as np
import h5py
from numpy import random
import sys

import astropy.units as u
import astropy.constants as const
import astropy.constants as c

import scipy.integrate
import scipy.interpolate

FloatType = np.float32
IntType = np.int32

class CentralStar(object):
    def __init__(self, M_star, Q, n, X_H, T_init, Z, Zi):
        self.M_star = M_star
        self.Q = Q
        self.n = n
        self.X_H = X_H
        self.n_H = n*X_H*4/(1+3*X_H)
        self.T_init = T_init
        self.Z = Z
        self.Zi = Zi
        self.aB = self.alpha_B()
        self.Rs = self.StromgrenRadius()
        self.trec = self.RecTime()
        self.StromgrenMass = self.IonizedMass()/self.X_H
        
    def alpha_B(self):
        return 2.59e-13*(1e4*u.K/(1e4*u.K))**(-0.833-0.034*np.log((self.T_init/(1e4*u.K))))*u.cm**3/u.s
    
    def StromgrenRadius(self):
        return (((3*self.Q)/(4*np.pi*(self.n_H)**2*self.aB))**(1/3)).to(u.pc)
        
    def RecTime(self):
        return (1/self.n_H/self.aB).to(u.Myr)   

    def IonizedMass(self):
        return (c.m_p/self.aB*self.Q/self.n_H).to(u.solMass)

class StromgrenSphere(object):
    def __init__(self, M_star, Q, n, X_H, T_init):
        self.M_star = M_star
        self.Q = Q
        self.n = n
        self.X_H = X_H
        self.n_H = n*X_H*4/(1+3*X_H)
        self.T_init = T_init
        self.aB = self.alpha_B()
        self.Rs = self.StromgrenRadius()
        self.trec = self.RecTime()
        self.StromgrenMass = self.IonizedMass()
        
    def alpha_B(self):
        return 2.59e-13*(1e4*u.K/(1e4*u.K))**(-0.833-0.034*np.log((self.T_init/(1e4*u.K))))*u.cm**3/u.s

    
    def StromgrenRadius(self):
        return (((3*self.Q)/(4*np.pi*(self.n_H)**2*self.aB))**(1/3)).to(u.pc)
        
    def RecTime(self):
        return (1/self.n_H/self.aB).to(u.Myr)   

    def IonizedMass(self):
        return (c.m_p/self.aB*self.Q/self.n_H).to(u.solMass)


def Write_to_hdf5(w, MassResolution, M_star, Z, Zi, BoxSize, CellsPerDimension, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal):
    #Type 4 particle
    Part4CellsPerDimension = 0
    Part4NumberOfCells = Part4CellsPerDimension**3+1
    #Part4MassPerCell = 1e-8
    Part4MassPerCell = M_star.value
    Part4Pos = np.zeros([Part4NumberOfCells, 3], dtype = FloatType)
    Random4Pos = random.rand(Part4NumberOfCells, 3)*BoxSize
    Part4Pos[:] = Random4Pos[:]
    Part4Mass = np.full(Part4NumberOfCells, Part4MassPerCell, dtype=FloatType)
    Part4Velocity = np.zeros([Part4NumberOfCells,3], dtype=FloatType)
    MassiveStarMass = np.full(Part4NumberOfCells, -1., dtype=FloatType)
    GFM_InitialMass = np.full(Part4NumberOfCells, Part4MassPerCell, dtype=FloatType)
    GFM_Metallicity = np.full((Part4NumberOfCells,),Z, dtype=FloatType)
    GFM_Metals = np.array(Zi, dtype=FloatType)
    GFM_StellarFormationTime = np.full((Part4NumberOfCells,), 1e-10, dtype=FloatType)
    Part4Pos[0,:] = np.array([BoxSize/2,BoxSize/2,BoxSize/2])
    MassiveStarMass[0] = M_star.value
    '''
    Part4CellsPerDimension = 12 
    Part4NumberOfCells = Part4CellsPerDimension**3+1 #1728 fake stars and one real star
    Part4MassPerCell = 1e-8
    Part4Pos = np.zeros([Part4NumberOfCells, 3], dtype = FloatType)
    Random4Pos = random.rand(Part4CellsPerDimension**3, 3)*BoxSize
    Part4Pos[1:,:] = Random4Pos[:,:]
    Part4Mass = np.full(Part4NumberOfCells, Part4MassPerCell, dtype=FloatType)
    Part4Velocity = np.zeros([Part4NumberOfCells,3], dtype=FloatType)
    MassiveStarMass = np.full(Part4NumberOfCells, 2, dtype=FloatType) #MassiveStarMass = 2 for fake stars
    GFM_InitialMass = np.zeros((Part4NumberOfCells,), dtype=FloatType)
    GFM_Metallicity = np.zeros((Part4NumberOfCells,), dtype=FloatType)
    GFM_Metals = np.zeros((Part4NumberOfCells, 9), dtype=FloatType)
    GFM_StellarFormationTime = np.full((Part4NumberOfCells,), 1e-10, dtype=FloatType)
    Part4Pos[0,:] = np.array([BoxSize/2,BoxSize/2,BoxSize/2])
    MassiveStarMass[0] = M_star.value
   '''
    if w == 0:
        if CellsPerDimension == 0:
            if bool(FIXED_RATE) == True:
                filename = "ICuB"+str(int(BoxSize))+'R'+str(int(MassResolution))+'Q'+str(int(FIXED_RATE))+'T'+str(int(InitialTemperature.value))+'.hdf5'
            if bool(FIXED_RATE) == False:
                filename = "ICuB"+str(int(BoxSize))+'R'+str(int(MassResolution))+'M'+str(int(M_star.value))+'T'+str(int(InitialTemperature.value))+'.hdf5'
        else:
            if bool(FIXED_RATE) == True:
                filename = "ICuB"+str(int(BoxSize))+'r'+str(int(CellsPerDimension))+'Q'+str(int(FIXED_RATE))+'T'+str(int(InitialTemperature.value))+'.hdf5'
            if bool(FIXED_RATE) == False:
                filename = "ICuB"+str(int(BoxSize))+'r'+str(int(CellsPerDimension))+'M'+str(int(M_star.value))+'T'+str(int(InitialTemperature.value))+'.hdf5'
    else :
        if bool(FIXED_RATE) == True:
            filename = "ICw"+str(w)+"B"+str(int(BoxSize))+'R'+str(int(MassResolution))+'Q'+str(int(FIXED_RATE))+'T'+str(int(InitialTemperature.value))+'.hdf5'
        if bool(FIXED_RATE) == False:
            filename = "ICw"+str(w)+"B"+str(int(BoxSize))+'R'+str(int(MassResolution))+'M'+str(int(M_star.value))+'T'+str(int(InitialTemperature.value))+'.hdf5'
    
    IC = h5py.File(filename, 'w')
    ## create hdf5 groups
    header = IC.create_group("Header")
    part0 = IC.create_group("PartType0")
    StarGroup = IC.create_group("PartType4")
    ## header entries
    NumPart = np.array([NumberOfCells, 0, 0, 0, Part4NumberOfCells, 0], dtype=IntType)
    header.attrs.create("BoxSize", BoxSize)
    header.attrs.create("Composition_vector_length", 0)
    header.attrs.create("Flag_Cooling", 0)
    header.attrs.create("Flag_DoublePrecision", 0)
    header.attrs.create("Flag_Feedback", 1)
    header.attrs.create("Flag_Metals", 0)
    header.attrs.create("Flag_Sfr", 1)
    header.attrs.create("Flag_StellarAge", 0)
    header.attrs.create("HubbleParam", 1.0)
    header.attrs.create("MassTable", np.zeros(6, dtype=IntType))
    header.attrs.create("NumFilesPerSnapshot", 1)
    header.attrs.create("NumPart_ThisFile", NumPart)
    header.attrs.create("NumPart_Total", NumPart)
    header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType))
    header.attrs.create("Omega0", 0.0)
    header.attrs.create("OmegaBaryon", 0.0)
    header.attrs.create("OmegaLambda", 0.0)
    header.attrs.create("Redshift", 0.0)
    header.attrs.create("Time", 0.0)
    header.attrs.create("UnitLength_in_cm", 3.085678e+18)
    header.attrs.create("UnitMass_in_g", 1.989e+33)
    header.attrs.create("UnitVelocity_in_cm_per_s", 10000.0)

    ## copy datasets
    part0.create_dataset("ParticleIDs", data=np.uint32(np.arange(1, NumberOfCells+1)))
    part0.create_dataset("Coordinates", data=Pos)
    part0.create_dataset("Masses", data=Mass)
    part0.create_dataset("Velocities", data=Velocity)
    part0.create_dataset("InternalEnergy", data=Uthermal)

    StarGroup.create_dataset("ParticleIDs", data=np.uint32(np.arange(NumberOfCells+1, 2*NumberOfCells+1)))
    StarGroup.create_dataset("Coordinates", data=Part4Pos)
    StarGroup.create_dataset("Velocities", data=Part4Velocity)
    StarGroup.create_dataset("Masses", data=Part4Mass)
    StarGroup.create_dataset("MassiveStarMass", data=MassiveStarMass)
    StarGroup.create_dataset("GFM_InitialMass", data=GFM_InitialMass)
    StarGroup.create_dataset("GFM_Metallicity", data=GFM_Metallicity)
    StarGroup.create_dataset("GFM_Metals", data=GFM_Metals)
    StarGroup.create_dataset("GFM_StellarFormationTime", data=GFM_StellarFormationTime)
    MassPerCell = Mass[0]
    ## close file
    IC.close()
    return filename

def Uniform_IC_Generate(centralStar, BoxSize, CellsPerDimension, FIXED_RATE, gamma = 5./3.):
    BoxSize = BoxSize.value
    NumberOfCells = int(CellsPerDimension**3)
    InitialDensity = centralStar.n.to(u.pc**(-3))
    InitialTemperature = centralStar.T_init
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    CellDensity = (InitialDensity*mu*c.m_p).to(u.solMass*u.pc**(-3)).value
    VolumePerCell = BoxSize**3/NumberOfCells
    MassPerCell = VolumePerCell*CellDensity
    DensityPerCell = MassPerCell/VolumePerCell
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value
    #Type 1 particle
    ID = np.zeros(NumberOfCells, dtype = IntType)
    Pos = np.zeros([NumberOfCells, 3], dtype = FloatType)
    Velocity = np.zeros([NumberOfCells, 3], dtype = FloatType)
    Mass = np.zeros(NumberOfCells, dtype = FloatType)
    Uthermal = np.zeros(NumberOfCells, dtype = FloatType)
    PhotonDensity = np.zeros(NumberOfCells, dtype = FloatType)
    RandomPos = random.rand(NumberOfCells, 3)*BoxSize
    Pos = RandomPos
    ## mass insetad of density
    Mass = np.full(NumberOfCells, MassPerCell, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)

    filename = Write_to_hdf5(0, 0, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, CellsPerDimension, NumberOfCells,  FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)

def Staggered_IC_Generate(centralStar, BoxSize, MassResolution, FIXED_RATE, gamma = 5./3.):
    InitialDensity = centralStar.n.to(u.pc**(-3))
    InitialTemperature = centralStar.T_init
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    CellDensity = (InitialDensity*mu*c.m_p).to(u.solMass*u.pc**(-3)).value
    BoxSize = BoxSize.value
    Mi = centralStar.StromgrenMass
    BoxMass = (((BoxSize*u.pc)**3*centralStar.n*c.m_p).to(u.solMass)).value
    TargetMass = (Mi/MassResolution).value
    NumberOfCells = BoxMass/(Mi.value/MassResolution)
    CellsPerDimension = int((NumberOfCells/2)**(1/3)+1)
    NumberOfCells = 2*CellsPerDimension**3
    BoxMass = NumberOfCells*TargetMass
    BoxSize = (BoxMass/CellDensity)**(1/3)
    VolumePerCell = BoxSize**3/NumberOfCells
    MassPerCell = TargetMass
    DensityPerCell = MassPerCell/VolumePerCell
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value    #Type 1 particle
    ## spacing
    dx = BoxSize / FloatType(CellsPerDimension)
    ## position of first and last cell
    pos_first, pos_last = 0.55/2 * dx, BoxSize - (1-0.55/2) * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
    Pos[:NumberOfCells//2,0] = xx.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,1] = yy.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,2] = zz.reshape(NumberOfCells//2)
    pos_first, pos_last = (0.55/2+0.45) * dx, BoxSize - (1-0.55/2-0.45) * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos[NumberOfCells//2:,0] = xx.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,1] = yy.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,2] = zz.reshape(NumberOfCells//2)
    Mass = np.full(NumberOfCells, MassPerCell, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)
    filename = Write_to_hdf5(0, MassResolution, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, 0, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)
    return filename

def StaggeredRandomize_IC_Generate(centralStar, BoxSize, MassResolution, FIXED_RATE, gamma = 5./3.):
    InitialDensity = centralStar.n.to(u.pc**(-3))
    InitialTemperature = centralStar.T_init
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    CellDensity = (InitialDensity*mu*c.m_p).to(u.solMass*u.pc**(-3)).value
    BoxSize = BoxSize.value
    Mi = centralStar.StromgrenMass
    BoxMass = (((BoxSize*u.pc)**3*centralStar.n*c.m_p).to(u.solMass)).value
    TargetMass = (Mi/MassResolution).value
    NumberOfCells = BoxMass/(Mi.value/MassResolution)
    CellsPerDimension = int((NumberOfCells/2)**(1/3)+1)
    NumberOfCells = 2*CellsPerDimension**3
    BoxMass = NumberOfCells*TargetMass
    BoxSize = (BoxMass/CellDensity)**(1/3)
    VolumePerCell = BoxSize**3/NumberOfCells
    MassPerCell = TargetMass
    DensityPerCell = MassPerCell/VolumePerCell
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value    #Type 1 particle
    ## spacing
    dx = BoxSize / FloatType(CellsPerDimension)
    ## position of first and last cell
    #pos_first, pos_last = 0.55/2 * dx, BoxSize - (1-0.55/2) * dx
    pos_first, pos_last = 0.5 * dx, BoxSize - 0.5 * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
    Pos[:NumberOfCells//2,0] = xx.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,1] = yy.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,2] = zz.reshape(NumberOfCells//2)
    #pos_first, pos_last = (0.55/2+0.45) * dx, BoxSize - (1-0.55/2-0.45) * dx
    pos_first, pos_last = (0.5-0.45) * dx, BoxSize - (0.5+0.45) * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos[NumberOfCells//2:,0] = xx.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,1] = yy.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,2] = zz.reshape(NumberOfCells//2)

    Delta_x = np.random.randn(Pos.shape[0],Pos.shape[1])*0.2*dx
    Pos = Pos + Delta_x
    Mass = np.full(NumberOfCells, MassPerCell, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)
    filename = Write_to_hdf5(0, MassResolution, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, 0, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)
    return filename

def Staggered_IC_TargetMass(centralStar, BoxSize, TargetMass, FIXED_RATE, gamma = 5./3.):
    InitialDensity = centralStar.n.to(u.pc**(-3))
    InitialTemperature = centralStar.T_init
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    CellDensity = (InitialDensity*mu*c.m_p).to(u.solMass*u.pc**(-3)).value
    TargetMass = (TargetMass.to(u.solMass)).value
    BoxSize = BoxSize.value
    BoxMass = (((BoxSize*u.pc)**3*centralStar.n*c.m_p).to(u.solMass)).value
    NumberOfCells = BoxMass/(TargetMass)
    CellsPerDimension = int((NumberOfCells/2)**(1/3)+1)
    NumberOfCells = 2*CellsPerDimension**3
    BoxMass = NumberOfCells*TargetMass
    BoxSize = (BoxMass/CellDensity)**(1/3)
    VolumePerCell = BoxSize**3/NumberOfCells
    MassPerCell = TargetMass
    DensityPerCell = MassPerCell/VolumePerCell
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value    #Type 1 particle
    ## spacing
    dx = BoxSize / FloatType(CellsPerDimension)
    ## position of first and last cell
    #pos_first, pos_last = 0.55/2 * dx, BoxSize - (1-0.55/2) * dx
    pos_first, pos_last = 0.5 * dx, BoxSize - 0.5 * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
    Pos[:NumberOfCells//2,0] = xx.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,1] = yy.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,2] = zz.reshape(NumberOfCells//2)
    #pos_first, pos_last = (0.55/2+0.45) * dx, BoxSize - (1-0.55/2-0.45) * dx
    pos_first, pos_last = (0.5-0.45) * dx, BoxSize - (0.5+0.45) * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos[NumberOfCells//2:,0] = xx.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,1] = yy.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,2] = zz.reshape(NumberOfCells//2)

    Delta_x = np.random.randn(Pos.shape[0],Pos.shape[1])*0.2*dx
    Pos = Pos + Delta_x
    Mass = np.full(NumberOfCells, MassPerCell, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)
    filename = Write_to_hdf5(0, TargetMass, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, 0, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)
    return filename

def Staggered_IC_Generate_Wind(centralStar, BoxSize, CellsPerDimension, FIXED_RATE, gamma = 5./3.):
    InitialDensity = centralStar.n.to(u.pc**(-3))
    InitialTemperature = centralStar.T_init
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    CellDensity = (InitialDensity*mu*c.m_p).to(u.solMass*u.pc**(-3)).value
    BoxSize = BoxSize.value
    Mi = centralStar.StromgrenMass
    BoxMass = (((BoxSize*u.pc)**3*centralStar.n*c.m_p).to(u.solMass)).value
    NumberOfCells = 2*CellsPerDimension**3
    VolumePerCell = BoxSize**3/NumberOfCells
    MassPerCell = VolumePerCell*CellDensity
    DensityPerCell = MassPerCell/VolumePerCell
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value    #Type 1 particle
    ## spacing
    dx = BoxSize / FloatType(CellsPerDimension)
    ## position of first and last cell
    pos_first, pos_last = 0.55/2 * dx, BoxSize - (1-0.55/2) * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
    Pos[:NumberOfCells//2,0] = xx.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,1] = yy.reshape(NumberOfCells//2)
    Pos[:NumberOfCells//2,2] = zz.reshape(NumberOfCells//2)
    pos_first, pos_last = (0.55/2+0.45) * dx, BoxSize - (1-0.55/2-0.45) * dx
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    Pos[NumberOfCells//2:,0] = xx.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,1] = yy.reshape(NumberOfCells//2)
    Pos[NumberOfCells//2:,2] = zz.reshape(NumberOfCells//2)
    Mass = np.full(NumberOfCells, MassPerCell, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)
    filename = Write_to_hdf5(0, CellsPerDimension, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, 0, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)
    return filename


def IC_Gen_Profile(profile, mPart):
    rProf = profile[0]
    rMax = np.max(rProf)
    vSphere = 4.*np.pi*rMax*rMax*rMax/3.
    rhoProf = profile[1]
    Mtot = scipy.integrate.simps(rhoProf*4.*np.pi*rProf*rProf, rProf)
    nPart = Mtot/mPart
    nProf = rhoProf/mPart
    n0 = nPart/vSphere

    # now only consider linear spacing case with the same dr
    dr = np.mean(rProf[1:]-rProf[:-1])

    # cumsum or integration
    rNew = ( np.cumsum(nProf*rProf*rProf*dr)/(n0/3.) ) ** (1/3.)
    rNew[-1] = rMax
    
    r_factor_inter = scipy.interpolate.interp1d( rNew, rProf/rNew )
    # uniform IC within a rMax**3 cube, sample particle number needs to be larger than nPart by Vcube/Vsphere
    samplePointsPerDim = int((nPart*8*rMax**3/vSphere)**(1/3))
    coordSample = np.zeros((samplePointsPerDim**3,3))
    Grid1d = np.linspace(-rMax, rMax, samplePointsPerDim, dtype=FloatType)
    xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
    coordSample[:,0] = xx.reshape(samplePointsPerDim**3)
    coordSample[:,1] = yy.reshape(samplePointsPerDim**3)
    coordSample[:,2] = zz.reshape(samplePointsPerDim**3)

    # random IC within a rMax**3 cube, sample particle number needs to be larger than nPart by Vcube/Vsphere
    #coordSample = np.random.uniform(-rMax,rMax,[int(nPart*8*rMax**3/vSphere),3])
    rSample = np.sqrt( np.sum( coordSample**2, axis=1 ) )
    coordSample = coordSample[rSample<rMax]
    rSample = rSample[rSample<rMax]
    
    rFactor = r_factor_inter(rSample)
    assert(rFactor.max()!=np.nan)
    
    coordSampleNew = coordSample*rFactor[:, None]
    
    return nPart, coordSampleNew

def Prof_ISO_core(r, rho0, r0, rMax, w):
    ind = (r < r0)
    ind0 = (r == 0)
    return rho0*ind + (1-ind)*rho0*((r+0.000001*ind0)/r0)**(-w)

def Prof_ISO_core_with_floor(r, rho0, r0, rMax, w, floorfactor):
    ind = (r < r0)
    ind0 = (r == 0)
    rholist = rho0*ind + (1-ind)*rho0*((r+0.000001*ind0)/r0)**(-w)
    rholist[rholist<floorfactor*rho0] = floorfactor*rho0
    return rholist

def Powerlaw_IC_Generate(centralStar, MassResolution, w, r0, BoxSize, FIXED_RATE, gamma = 5./3.):
    BoxSize = BoxSize.value
    r0 = r0.value
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    rho0 = ((mu*c.m_p*centralStar.n).to(u.solMass*u.pc**(-3))).value
    Mi = centralStar.StromgrenMass
    TargetMass = (Mi/MassResolution).value
    rMax = BoxSize/2*np.sqrt(3)
    rIsoCore = np.linspace(0,rMax,1000000)
    profIsoCore = Prof_ISO_core(rIsoCore, rho0, r0, rMax, w)
    nPart, coord = IC_Gen_Profile([rIsoCore, profIsoCore], TargetMass)
    coord+=BoxSize/2
    ind_list = []
    for i in range(coord.shape[0]):
        if (0<coord[i,0]<BoxSize) and (0<coord[i,1]<BoxSize) and (0<coord[i,2]<BoxSize):
            ind_list.append(i)
    coord = coord[ind_list]
    NumberOfCells = coord.shape[0]
    InitialTemperature = centralStar.T_init
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value    
    Pos = np.zeros([NumberOfCells, 3], dtype = FloatType)
    Velocity = np.zeros([NumberOfCells, 3], dtype = FloatType)
    Mass = np.zeros(NumberOfCells, dtype = FloatType)
    Uthermal = np.zeros(NumberOfCells, dtype = FloatType)
    Pos = coord
    ## mass insetad of density
    Mass = np.full(NumberOfCells, TargetMass, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)
    filename = Write_to_hdf5(w, MassResolution, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, 0, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)
    return filename

def Powerlawfloor_IC_Generate(centralStar, MassResolution, w, r0, BoxSize, floorfactor, FIXED_RATE, gamma = 5./3.):
    BoxSize = BoxSize.value
    r0 = r0.value
    X_H = centralStar.X_H 
    mu = 4/(1+3*X_H)
    rho0 = ((mu*c.m_p*centralStar.n).to(u.solMass*u.pc**(-3))).value
    Mi = centralStar.StromgrenMass
    TargetMass = (Mi/MassResolution).value
    rMax = BoxSize/2*np.sqrt(3)
    rIsoCore = np.linspace(0,rMax,1000000)
    profIsoCore = Prof_ISO_core_with_floor(rIsoCore, rho0, r0, rMax, w, floorfactor)
    nPart, coord = IC_Gen_Profile([rIsoCore, profIsoCore], TargetMass)
    coord+=BoxSize/2
    ind_list = []
    for i in range(coord.shape[0]):
        if (0<coord[i,0]<BoxSize) and (0<coord[i,1]<BoxSize) and (0<coord[i,2]<BoxSize):
            ind_list.append(i)
    coord = coord[ind_list]
    NumberOfCells = coord.shape[0]
    InitialTemperature = centralStar.T_init
    UthermalPerCell = (c.k_B*InitialTemperature/mu/c.m_p/(gamma-1)).to((u.km/u.s)**2).value    
    Pos = np.zeros([NumberOfCells, 3], dtype = FloatType)
    Velocity = np.zeros([NumberOfCells, 3], dtype = FloatType)
    Mass = np.zeros(NumberOfCells, dtype = FloatType)
    Uthermal = np.zeros(NumberOfCells, dtype = FloatType)
    Pos = coord
    ## mass insetad of density
    Mass = np.full(NumberOfCells, TargetMass, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Uthermal = np.full(NumberOfCells, UthermalPerCell, dtype=FloatType)
    filename = Write_to_hdf5(w, MassResolution, centralStar.M_star, centralStar.Z, centralStar.Zi, BoxSize, 0, NumberOfCells, FIXED_RATE, InitialTemperature, Pos, Mass, Velocity, Uthermal)
    return filename

