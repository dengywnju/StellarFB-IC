import numpy as np
import h5py
from numpy import random
import Stromgren_IC_utils as IC
import astropy.units as u
import astropy.constants as const
import astropy.constants as c
import IC_Checker

#######Select a type of IC###########################################################
########Mass Resolution##############################################################
#TypeOfIC = 'Uniform-Regular'          #uniform density medium with regular Cartesian mesh
#TypeOfIC = 'Uniform-Staggered'         #uniform density medium with staggered Cartesian mesh
TypeOfIC = 'Uniform-TargetMass'
#TypeOfIC = 'Uniform-Staggered-randomize'         #uniform density medium with staggered Cartesian mesh with 0.2Dx gaussian offset
#TypeOfIC = 'Powerlaw'                 #uniform core + powlaw skirt
#TypeOfIC = 'Powerlaw_with_floor'      #uniform core + powlaw skirt + uniform background
#######Spatial Resolution#############################################################
#TypeOfIC = 'Uniform-Staggered-wind'   #uniform density medium with staggered Cartesian mesh

#####Set Initial Conditions###########################################################
#Star
SolarMetallicity = True

M_Star = 20*u.solMass #initial stellar mass
Q = 1e48/u.s #total ionising photon rate (for FIXED_RATE mode)
Q0 = (1e63/u.Gyr).to(1/u.s)
FIXED_RATE = (Q/Q0).value

#Gas
BoxSize = 20*u.pc
gamma = 5./3.
n = 100/u.cm**3  #core density for the powerlaw case
T_init = 1e4*u.K
if SolarMetallicity == True:
    Z = 0.0134 #Metallicity
    Zi = np.array([0.7381,0.2485,2.4e-3,7e-4,5.8e-3,1.3e-3,7e-4,7e-4,1.3e-3,5e-4]) #metals C, N, O, Ne, Mg, Si, S, Ca, Fe, and other metals.
    X_H = Zi[0]
else:
    X_H = 1 #H abundance, 0.76 for primordial, 1 for pure hydrogen
    Z = 1e-2 #Metallicity
    Zi = np.array([5e-3,3e-3,2e-3,0,0,0,0,0,0]) #metals C, N, O, Ne, Mg, Si, S, Ca, Fe


#Resolution
CellsPerDimension = 64 #Spatial resolution, !! for Uniform-Staggered-wind mode, total number of cells is 2*CellsPerDimension**3
MassResolution = 10000 #Ionized Mass resolution
TargetMass = 1*u.solMass #TargetMass mode

#Specific parameters for powerlaw modes
w = 2
r0 = 0.75*u.pc
floorfactor = 1e-4
####################################################################################
#stromgren = IC.StromgrenSphere(M_Star, Q, n, X_H, T_init)
centralStar = IC.CentralStar(M_Star, Q, n, X_H, T_init, Z, Zi)


if TypeOfIC == 'Uniform-Regular':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Gas Density:", str(n.value)+' /cm3')
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: CellsPerDimension:", CellsPerDimension)
    filename = IC.Uniform_IC_Generate(centralStar, BoxSize, 
         CellsPerDimension, int(FIXED_RATE), gamma)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

if TypeOfIC == 'Uniform-Staggered-randomize':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Gas Density:", str(n.value)+' /cm3')
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: MassResolution:", MassResolution)
    filename = IC.StaggeredRandomize_IC_Generate(centralStar, BoxSize, 
         MassResolution, int(FIXED_RATE), gamma)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

if TypeOfIC == 'Uniform-Staggered':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Gas Density:", str(n.value)+' /cm3')
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: MassResolution:", MassResolution)
    filename = IC.Staggered_IC_Generate(centralStar, BoxSize, 
         MassResolution, int(FIXED_RATE), gamma)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

if TypeOfIC == 'Uniform-Staggered-wind':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Gas Density:", str(n.value)+' /cm3')
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: CellsPerDimension:", CellsPerDimension)
    filename = IC.Staggered_IC_Generate_Wind(centralStar, BoxSize, 
         CellsPerDimension, int(FIXED_RATE), gamma)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

if TypeOfIC == 'Uniform-TargetMass':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Gas Density:", str(n.value)+' /cm3')
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: CellsPerDimension:", CellsPerDimension)
    filename = IC.Staggered_IC_TargetMass(centralStar, BoxSize, 
         TargetMass, int(FIXED_RATE), gamma)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)


if TypeOfIC == 'Powerlaw':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: Radius of Core:", str(r0.value)+' pc')
    print("IC-Generator: Core Density:", str(n.value)+' /cm3')
    print("IC-Generator: Power-law index:", w)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: CellsPerDimension:", MassResolution)
    filename = IC.Powerlaw_IC_Generate(centralStar, MassResolution, w, r0, BoxSize, 
     int(FIXED_RATE), gamma = 5./3.)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

if TypeOfIC == 'Powerlaw_with_floor':
    print("IC-Generator: Generating...")
    print("IC-Generator: Type of IC:", TypeOfIC)
    print("IC-Generator: Floor factor:", str(floorfactor))
    print("IC-Generator: Gas Temperature:", str(T_init.value)+' K')
    print("IC-Generator: gamma:", gamma)
    print("IC-Generator: Radius of Core:", str(r0.value)+' pc')
    print("IC-Generator: Core Density:", str(n.value)+' /cm3')
    print("IC-Generator: Power-law index:", w)
    print("IC-Generator: MassiveStarMass:", str(M_Star.value)+' Msun')
    print("IC-Generator: FIXED_RATE", FIXED_RATE)
    print("IC-Generator: BoxSize:", str(BoxSize.value)+' pc')
    print("IC-Generator: MassResolution:", MassResolution)
    filename = IC.Powerlawfloor_IC_Generate(centralStar, MassResolution, w, r0, BoxSize, floorfactor, int(FIXED_RATE), gamma = 5./3.)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

