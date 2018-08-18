#!/usr/bin/env python
from __future__ import print_function

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt

import PlotScripts
import ReadScripts
import AllVars

alpha = 6  # Best fit from Kravstov, Gnedin, Klypin 2004.

def get_jeans_mass(hubble_h, omega_m):

    mu = 0.59  # Mean molecular weight of fully ionized gas. 
    jeans_mass = 2.5e11 / hubble_h * pow(omega_m, -0.5) * pow(mu, -1.5) 

    return jeans_mass


def get_analytic_fit(Redshift, z_0, z_r):

    a_0 = 1.0 / (1.0 + z_0)  # Analytic function uses redshift. 
    a_r = 1.0 / (1.0 + z_r)
    a = 1.0 / (1.0 + Redshift)
    try:
        f = np.empty(len(a))
    except TypeError:
        is_array = False 
    else:
        is_array = True 

    if is_array:
        for count, scale_factor in enumerate(a):
            if scale_factor < a_0:
                f[count] = 3.0 * scale_factor / ((2+alpha)*(5+2*alpha)) \
                           * pow(scale_factor / a_0, alpha)
            elif scale_factor > a_0 and scale_factor < a_r:
                f[count] = 3.0 / scale_factor * (a_0*a_0* (1 / (2+alpha) - 2*pow(scale_factor / a_0, -0.5) / (5 + 2*alpha)) \
                                                +scale_factor*scale_factor / 10.0 \
                                                -a_0*a_0 / 10.0 * (5 - 4*pow(scale_factor / a_0, -0.5))) 
                                                 
            else:
                f[count] = 3.0 / scale_factor * (a_0*a_0* (1 / (2+alpha) - 2*pow(scale_factor / a_0, -0.5) / (5 + 2*alpha)) \
                                                +a_r*a_r / 10.0 * (5 - 4*pow(scale_factor / a_0, -0.5)) \
                                                -a_0*a_0 / 10.0 * (5 - 4*pow(scale_factor / a_0, -0.5)) \
                                                +scale_factor*a_r / 3.0 \
                                                -a_r*a_r / 3.0 * (3 - 2*pow(scale_factor / a_r, -0.5))) 

    else:
        scale_factor = a 
        if scale_factor < a_0:
            f = 3.0 * scale_factor / ((2+alpha)*(5+2*alpha)) \
                       * pow(scale_factor / a_0, alpha)
        elif scale_factor > a_0 and scale_factor < a_r:
            f = 3.0 / scale_factor * (a_0*a_0* (1 / (2+alpha) - 2*pow(scale_factor / a_0, -0.5) / (5 + 2*alpha)) \
                                            +scale_factor*scale_factor / 10.0 \
                                            -a_0*a_0 / 10.0 * (5 - 4*pow(scale_factor / a_0, -0.5))) 
                                             
        else:
            f = 3.0 / scale_factor * (a_0*a_0* (1 / (2+alpha) - 2*pow(scale_factor / a_0, -0.5) / (5 + 2*alpha)) \
                                            +a_r*a_r / 10.0 * (5 - 4*pow(scale_factor / a_0, -0.5)) \
                                            -a_0*a_0 / 10.0 * (5 - 4*pow(scale_factor / a_0, -0.5)) \
                                            +scale_factor*a_r / 3.0 \
                                            -a_r*a_r / 3.0 * (3 - 2*pow(scale_factor / a_r, -0.5))) 

    return f

def get_filter_mass(Redshift, z_0, z_r):

    jeans_mass = get_jeans_mass(AllVars.Hubble_h, AllVars.Omega_m) 
    
    f = get_analytic_fit(Redshift, z_0, z_r)

    filter_mass = np.log10(jeans_mass * pow(f, 1.5))

    '''
    for count, z in enumerate(Redshift):
        print("For redshift {0:.2f} filter mass is {1:.4e}".format(z, filter_mass[count])) 
    '''

    return filter_mass


def plot_reionmod(HaloMass, FilterMass, Redshift, fb):

    plot_redshift = [10, 9, 8, 7, 6]
  
    fig1 = plt.figure() 
    ax1 = fig1.add_subplot(111)
    PlotScripts.Set_Params_Plot()
 
    for count, z in enumerate(plot_redshift):
        z_idx = np.abs(Redshift - z).argmin()
        
        reionmod = 1 / pow(1 + 0.26*pow(10, FilterMass[z_idx] - HaloMass), 3.0)
 
        label = "z = {0:d}".format(z)       
        ax1.plot(HaloMass, reionmod, color = PlotScripts.colors[count], 
                 label = label)

    ax1.set_xlim([7.0, 12.0])
    ax1.set_ylim([0.0, 1.05])

    ax1.set_xlabel(r'$\mathbf{log_{10} \: M_{vir} \:[M_{\odot}]}$', fontsize = PlotScripts.global_labelsize) 
    ax1.set_ylabel(r'$\mathbf{ReionMod}$', fontsize = PlotScripts.global_labelsize) 

    leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)

    outputFile1 = "./gnedin.png"
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


if __name__ == "__main__":

    HaloMass = np.arange(7, 12, 0.1)
    Redshift = np.arange(15, 6, -0.1)
    z_0 = 8
    z_r = 7
    AllVars.Set_Params_Kali()

    FilterMass = get_filter_mass(Redshift, z_0, z_r)
    plot_reionmod(HaloMass, FilterMass, Redshift, AllVars.BaryonFrac)
