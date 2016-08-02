from scipy import optimize
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

R = 1.986e-3
CtoK = 273.15
def melt_curve(T, Tm, dH, min, max):
    """ functional form of a melting curve, parametrized by Tm and dH """
    K = np.exp((dH/R) * (1/(Tm+CtoK) - 1/(T+CtoK)))
    return min + max*1/(1+K)

def fit_melt_curve(T, p, sigma):
    """
    fit melting curve to given array of temperatures and measured proportion
    folded
    returns fit parameters and covariance matrix
    """
    return optimize.curve_fit(melt_curve, T, p, p0=[60, -35, 0, 1],
                              sigma=sigma)

def melt_curves(T, Tm, dH, *args):
    """
    functional form of a melting curve, with multiple min/max choices
    inputs:
        x = [T, choice] where T is the temperature and
                              choice is the choice of min/max
        Tm = melting temperature
        dH = enthalpy
        args should contain minX and maxX for all choices X
    """
    K = np.exp((dH/R) * (1/(T[0,:]+CtoK) - 1/(Tm+CtoK)))
    min = np.array([args[int(x)*2] for x in T[1,:]])
    max = np.array([args[int(x)*2+1] for x in T[1,:]])
    return min + max*1/(1+K)

def fit_melt_curves(T, p, sigma):
    """
    fit melting curve to given array of temperatures and measured proportion
    folded
    returns fit parameters and covariance matrix
    """
    ncurves = p.shape[1]
    p = p.values.flatten('F') 
    sigma = sigma.values.flatten('F')
    newT = []
    for i in range(ncurves):
        newT.extend([[t, i] for t in T])
    newT = np.array(newT).T
    p0 = [60, -35] + [0, 1]*ncurves
    return optimize.curve_fit(melt_curves, newT, p, p0=p0, sigma=sigma)

def fit_and_plot_melts(T, p, sigma):
    """ fit melting curves and plot data """
    popt, pcov = fit_melt_curves(T, p, sigma)
    for i in range(p.shape[1]):
        p.ix[:,i] -= popt[(i+1)*2]
        p.ix[:,i] /= popt[(i+1)*2+1]
        sigma.ix[:,i] /= popt[(i+1)*2+1]
        plt.errorbar(T, 1-p.ix[:,i], sigma.ix[:,i], fmt='o',
                     label=p.columns[i])
    plot_melt_curve([popt[0], popt[1], 0, 1])
    plt.text(20, plt.gca().get_ylim()[1]*.9,
             r'$T_m = %.1f\pm %.1f, \Delta H = %.1f\pm %.1f$' %
             (popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])))
    plt.xlabel(r'temperature ($^\circ$C)')
    plt.ylabel('reactivity')
    plt.legend()
    return popt, pcov

def plot_melt_curve(params, Ts=range(100)):
    """ plot melt curve with given params over given temperatures """
    p = [melt_curve(T, params[0], params[1], params[2], params[3]) for T in Ts]
    plt.plot(Ts, p, 'k-')
    return

def fit_and_plot_melt(T, p, sigma):
    """ fit and plot data """
    popt, pcov = fit_melt_curve(T, p, sigma)
    plot_melt_curve(popt)
    plt.errorbar(T, p, sigma, fmt='o')
    plt.text(20, plt.gca().get_ylim()[1]*.9,
             r'$T_m = %.1f\pm %.1f, \Delta H = %.1f\pm %.1f$' %
             (popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])))
    plt.xlabel(r'temperature ($^\circ$C)')
    plt.ylabel('reactivity')
    return popt, pcov
        
def parse_rdat(filename):
    """ parse data from rdat file """
    Ts = []
    cols = []
    data = []
    err = []
    with open(filename) as f:
        for line in f:
            if line.startswith('ANNOTATION_DATA'):
                Ts.append(float(line.split('\t')[1].strip('temperature: ')))
            elif line.startswith('SEQPOS'):
                cols = line.split()[1:]
            elif line.startswith('REACTIVITY:'):
                data.append([float(x) for x in line.split()[1:]]) 
            elif line.startswith('REACTIVITY_ERROR:'):
                err.append([float(x) for x in line.split()[1:]]) 
    data = pd.DataFrame(data, columns=cols)
    err = pd.DataFrame(err, columns=cols)
    return Ts, data, err

def scale(data, err):
    """ scale data to range [0,1] and error accordingly """
    data -= data.min()
    max = data.max()
    data /= max
    err /= max
    return data, err

def normalize(data, err, bases):
    """ normalize data to given bases """
    data, err = scale(data, err)
    norm = data.ix[:,bases].mean(axis=1)
    data = data.divide(norm, axis='index')
    err = err.divide(norm, axis='index')
    return data    
