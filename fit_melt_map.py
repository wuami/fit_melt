import os
import fit_utils as u
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

results_dir = RESULTSDIR 

def fit_bases(T, data, err, bases, title, suffix=''):
    fits = pd.DataFrame(index=['%s%s' % (base, suffix) for base in bases],
                        columns=['Tm', 'dH', 'Tm_stderr', 'dH_stderr'])
    for col in bases:
        try:
            plt.figure()
            popt, pcov = u.fit_and_plot_melt(T, data.ix[:,col], err.ix[:,col])
            plt.savefig('%s/%s_%s%s.png' % (results_dir, title, col, suffix))
            plt.close()
            fits.ix['%s%s' % (col, suffix),:] = [popt[0], popt[1], np.sqrt(pcov[0,0]),
                              np.sqrt(pcov[1,1])]
        except RuntimeError as e:
            plt.figure()
            plt.errorbar(T, data.ix[:,col], err.ix[:,col], fmt='o')
            plt.savefig('%s/%s_%s%s.png' % (results_dir, title, col, suffix))
            plt.close()
    return fits

def jointly_fit_bases(T, data, err, bases, title, suffix=''):
    name = '%s%s' % ('.'.join(bases), suffix)
    plt.figure()
    popt, pcov = u.fit_and_plot_melts(T, data.ix[:,bases], err.ix[:,bases])
    print name, popt
    plt.savefig('%s/%s_%s.png' % (results_dir, title, name))
    plt.close()
    result = {'Tm': popt[0], 'dH': popt[1], 'Tm_stderr': np.sqrt(pcov[0,0]),
              'dH_stderr': np.sqrt(pcov[1,1])}
    return pd.DataFrame(result, index=[name])

def parse_args():
    p = argparse.ArgumentParser(description='fit a melt curve')
    p.add_argument('-f', '--filename', help='name of the file containing data')
    p.add_argument('-o', '--outfile', help='name of output file', default=None)
    p.add_argument('-t', '--title', default=None,
                   help='title for plot and figure filenames')
    p.add_argument('-b', '--bases', nargs='+', help='positions to fit to')
    p.add_argument('-n', '--normalize', nargs='*', default=[], 
                   help='positions to normalize to')
    return p.parse_args()
                                        

def main():
    args = parse_args()
    T, data, err = u.parse_rdat(args.filename)

    # set defaults for optional arguments
    if args.title is None:
        args.title = os.path.basename(args.filename)
    if args.outfile is None:
        args.outfile = args.filename.rstrip('.rdat')
    
    # normalize data if necessary
    if len(args.normalize) > 0:
        data = u.normalize(data, err, args.normalize)

    # fit data
    if len(args.bases) == 1:
        fits = fit_bases(T, data, err, args.bases, args.title)
    elif len(args.bases) > 1:
        fits = jointly_fit_bases(T, data, err, args.bases, args.title)
    
    # write to file
    fits.to_csv('%s.fit' % args.outfile, sep='\t')


if __name__ == '__main__':
    main()
