import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../../')

import matplotlib.pyplot as plt
import scipy
from copy import deepcopy
from uncertainties import ufloat as uf
from uncertainties import unumpy as unp
from scipy import stats as spys
from itertools import chain
import numpy as np
import math

import HardRockMC.ExceltoPandasMethods as EPM
from fitter import Fitter, get_common_distributions, get_distributions

from HardRockMC.histfitter import histfitter
from HardRockMC.histplotter import histplotter


class PorosityScale:

    def allBH(filename,
      xparam='volume [ccm/piece]',
      xlabel=r'Volume of each core piece [cm$^{3}$]',
      expt_means=False,
      expt_legend=False):

        yparam = 'neighbor porosity'

        BHs = ['BH01', 'BH02', 'BH03']
        colors = ['tab:red', 'tab:olive', 'tab:purple']

        expts = ['PWoutdiff', 'AQ', 'IsoEx LAB', 'IsoEx ICE', 'BD', 'TP']
        markers = ['o', 'o', 'o', 'o', 'o', 'o', 'o']

        VolRanges = [(5e-1,7e-1), (3,10), (3e1,2e2), (3e2,1e3)]

        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(8,6))

        histplotter.regression_samplescale(axs, xparam, yparam,
          fitter_kwargs={'BH':'all', 'exclude_expts':['PP']},
          plot_kwargs={'scatter':False, 'label':'Regression'})

        x_scat = []
        y_scat = []
        for i, BH in enumerate(BHs):

            for expt, mk in zip(expts, markers):

                # this returns the x and y means and standards if mean is passed true.
                # otherwise, returns the x and y lists of data points
                x, xs, y, ys = histplotter.scatter_samplescale(
                  axs, xparam, yparam,
                  errorbars=True,
                  # mean=True,
                  fitter_kwargs={'BH':BH, 'only_expt':expt},
                  plot_kwargs={'marker':mk, 'color':colors[i],  'alpha':0.5,'zorder':0, 'ecolor':colors[i]})

                if expt_means:
                    axs.scatter(x,y, marker='*', s=10)

                if expt_legend:
                    if BH == BHs[-1]:
                        axs.scatter([np.nan], [np.nan], marker=mk, s=5., color='k', label=expt)

                # collect the xs and ys for the local volume mean and std dev
                [x_scat.append(_x) for _x in x]
                [y_scat.append(_y) for _y in y]

                # for BH color legend
                if expt == expts[0]:
                    axs.scatter([np.nan], [np.nan], marker="o", s=30., color=colors[i], label=BH)

            axs.set_ylabel('Primary porosity [vol/vol]')
            axs.set_xlabel(xlabel)
            axs.set_xscale('log')

        x_scat = np.array(x_scat)
        y_scat = np.array(y_scat)


        # add local mean for the volume ranges
        for V_tup in VolRanges:
            # the where sets anything higher than the upper limit to 0
            # then the asarray.nonzero selects the indicies that are higher than
            # the lower limit.
            # so these are the indicies of x and y where y is in this range.
            _indicies = np.asarray(np.where(x_scat < V_tup[1], x_scat, 0)>V_tup[0]).nonzero()

            this_x = x_scat[_indicies]
            this_y = y_scat[_indicies]

            bigeb = axs.errorbar(this_x.mean(), this_y.mean(),
              yerr=3*this_y.std(), marker='*', color='k', linewidth=2.0,
              capsize=6, markeredgewidth=2.)

            bigeb[-1][0].set_linestyle('--')

            axs.errorbar(this_x.mean(), this_y.mean(),
              yerr=this_y.std(), marker='*', markersize=15, color='k',
              linewidth=2.0, capsize=6, markeredgewidth=2.)

            axs.text(0.5*this_x.mean(), 0.00125, 'n = '+str(len(this_x))+'\n'+r'$\mu$ = '+str(round(this_y.mean(),4))+'\n'+r'$\sigma$ = '+str(round(this_y.std(),4)), va='top', ha='left')


        axs.scatter([np.nan], [np.nan], marker="*", color='k', label='Local mean')

        plt.legend(loc='upper right')

        plt.ylim(0,0.01)
        plt.xlim(2e-1,1e3)

        fig.suptitle('Primary porosity (based on core water content)')

        plt.tight_layout()
        plt.savefig(os.path.dirname(__file__)+'/../../Figures/'+filename)
        plt.close()
