import numpy as np
import h5py
from numpy import random
import Stromgren_IC_utils as IC
import astropy.units as u
import astropy.constants as const
import astropy.constants as c
import IC_Checker

#######Select a type of IC###########################################################

#TypeOfIC = 'Uniform-Regular'
#TypeOfIC = 'Uniform-Staggered'
#TypeOfIC = 'Powerlaw'
#TypeOfIC = 'Powerlaw_with_floor'
TypeOfIC = 'Uniform-Staggered-wind'

#####Set Initial Conditions###########################################################
#Star
M_Star = 4*u.solMass
Q = 1e48/u.s
Q0 = (1e63/u.Gyr).to(1/u.s)
FIXED_RATE = (Q/Q0).value

#Gas
BoxSize = 50*u.pc
gamma = 5./3.
n = 1/u.cm**3  #core density for the powerlaw case
T_init = 1e4*u.K
X_H = 0.76 #H abundance, 0.76 for primordial, 1 for pure hydrogen

#Resolution
CellsPerDimension = 64 #Spatial resolution, !! for Uniform-Staggered-wind mode, total number of cells is 2*CellsPerDimension**3
MassResolution = 1000 #Ionized Mass resolution

#Specific parameters for powerlaw modes
w = 2
r0 = 0.75*u.pc
floorfactor = 1e-4
####################################################################################
stromgren = IC.StromgrenSphere(M_Star, Q, n, X_H, T_init)

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
    filename = IC.Uniform_IC_Generate(stromgren, BoxSize, 
         CellsPerDimension, int(FIXED_RATE), gamma)
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
    filename = IC.Staggered_IC_Generate(stromgren, BoxSize, 
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
    filename = IC.Staggered_IC_Generate_Wind(stromgren, BoxSize, 
         CellsPerDimension, int(FIXED_RATE), gamma)
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
    filename = IC.Powerlaw_IC_Generate(stromgren, MassResolution, w, r0, BoxSize, 
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
    filename = IC.Powerlawfloor_IC_Generate(stromgren, MassResolution, w, r0, BoxSize, floorfactor, int(FIXED_RATE), gamma = 5./3.)
    print("IC-Generator: Filename:", filename)
    IC_Checker.Checker(filename, TypeOfIC, w, n, X_H, r0, gamma)

