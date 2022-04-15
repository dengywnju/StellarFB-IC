import matplotlib
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate
import scipy.interpolate
import astropy.units as u
import astropy.constants as c
import Stromgren_IC_utils as IC
plt.style.use('stylesheet.mplstyle')

def Checker(filename, ICType, w, n0, X_H, r0, gamma = 5./3):
    print('IC-Checker: Checking...')
    with h5py.File(filename,'r') as f:
        BoxSize = f['Header'].attrs['BoxSize']
        x_pos = f['PartType0']['Coordinates'][:,0]
        y_pos = f['PartType0']['Coordinates'][:,1]
        z_pos = f['PartType0']['Coordinates'][:,2]
        x_pos = x_pos-BoxSize/2
        y_pos = y_pos-BoxSize/2
        z_pos = z_pos-BoxSize/2
        uthermal = f['PartType0']['InternalEnergy'][:]
        mass = f['PartType0']['Masses'][:]

    Radius = np.sqrt(x_pos**2+y_pos**2+z_pos**2)
    MassPerCell = mass[0]
    NumberOfCells = mass.shape[0]
    CellsPerDimension = NumberOfCells**(1/3)
    mu = 4/(1+3*X_H)
    Temperature = ((mu*c.m_p*(gamma-1)*uthermal[0]*(u.km/u.s)**2/c.k_B).to(u.K)).value
    hist, edge = np.histogram(Radius, bins = int(CellsPerDimension), range = [0.05*BoxSize,BoxSize/2])
    rPlot = 0.5*(edge[:-1]+edge[1:])
    vShell = 4*np.pi/3.*(edge[1:]**3-edge[:-1]**3)
    Rho = hist*MassPerCell/vShell
    DensFac = (u.solMass/u.pc**3/(mu*c.m_p)).to(1/u.cm**3).value

    if ICType == 'Powerlaw':
        r0 = r0.value
        rho0 = ((n0*c.m_p).to(u.solMass/u.pc**3)).value
        rIsoCore = np.linspace(0,BoxSize/2,100000)
        profIsoCore = IC.Prof_ISO_core(rIsoCore, rho0, r0, 5*r0, w)
        profIsoCoreTest = IC.Prof_ISO_core(rPlot, rho0, r0, 5*r0, w)
        chi2 = np.sum((Rho-profIsoCoreTest)**2/profIsoCoreTest)

        print('IC-Checker: BoxSize:', BoxSize)
        print('IC-Checker: MassPerCell:', MassPerCell)
        print('IC-Checker: NumberOfCells:', NumberOfCells)
        print('IC-Checker: H abundance:', str(X_H))
        print('IC-Checker: Molecular Weight:', str(mu))
        print('IC-Checker: Temperature:', str(Temperature)+" K")
        print('IC-Checker: chi-square:', str(float('%.4g' % chi2)))
        print('IC-Checker: IC-Checker.pdf')

        plt.step(rPlot,Rho*DensFac, where = 'mid')
        plt.plot(rIsoCore, profIsoCore*DensFac, color = 'grey', linestyle ='--')
        plt.text(0.55,0.9, 'Core Density: '+ str(int(n0.value))+r" cm$^{-3}$",transform = plt.gca().transAxes)
        plt.text(0.55,0.85, 'Temperature: '+ str(int(Temperature))+" K",transform = plt.gca().transAxes)
        plt.text(0.55,0.8, 'MassPerCell: '+ str(float('%.6g' % MassPerCell)),transform = plt.gca().transAxes)
        plt.text(0.55,0.75, 'NumberOfCells: '+ str(format(NumberOfCells,'.2E')),transform = plt.gca().transAxes)
        plt.text(0.55,0.7, r'$\chi^2$: '+ str(float('%.4g' % chi2)),transform = plt.gca().transAxes)
        plt.xlim(rPlot[0], rPlot[-1])
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel(r'$R$')
        plt.ylabel(r'$\rho(R)$')  
        plt.savefig('IC-Checker.pdf')

    else:
        Density_ = np.sum(rPlot**2*Rho)/np.sum(rPlot**2)
        Density = Density_*DensFac

        print('IC-Checker: BoxSize:', BoxSize)
        print('IC-Checker: MassPerCell:', MassPerCell)
        print('IC-Checker: NumberOfCells:', NumberOfCells)
        print('IC-Checker: H abundance:', str(X_H))
        print('IC-Checker: Molecular Weight:', str(mu))
        print('IC-Checker: Temperature:', str(Temperature)+" K")
        print('IC-Checker: Density:', str(Density)+" cm-3")
        print('IC-Checker: IC-Checker.pdf')

        plt.step(rPlot,Rho*DensFac, where = 'mid')
        plt.axhline(Density_*DensFac, color = 'grey', linestyle ='--')
        plt.text(0.55,0.9, 'Density: '+ str(int(Density))+r" cm$^{-3}$",transform = plt.gca().transAxes)
        plt.text(0.55,0.85, 'Temperature: '+ str(int(Temperature))+" K",transform = plt.gca().transAxes)
        plt.text(0.55,0.8, 'MassPerCell: '+ str(float('%.6g' % MassPerCell)),transform = plt.gca().transAxes)
        plt.text(0.55,0.75, 'NumberOfCells: '+ str(format(NumberOfCells,'.2E')),transform = plt.gca().transAxes)
        plt.xlabel(r'$R$')
        plt.ylabel(r'$\rho(R)$')  
        plt.savefig('IC-Checker.pdf')

        
    
