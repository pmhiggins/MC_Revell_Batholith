import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')
sys.path.append(os.path.dirname(__file__)+'/../../')
import string


import ExceltoPandasMethods as EPM
import matplotlib.pyplot as plt
import scipy
from copy import deepcopy
from uncertainties import ufloat as uf
from uncertainties import unumpy as unp


from scipy import stats as spys
from itertools import chain
import numpy as np
import math
import pandas as pd

# import astropy.stats as apys
from fitter import Fitter, get_common_distributions, get_distributions

from histfitter import histfitter
from histplotter import histplotter

class BHbyBH_inputcomparison:

    xlimsdict = {
      'neighbor porosity' : [0.0, 0.014],
      'bulk density' : [2.5,3.2],
      'U' : [-0.5, 1.5],
      'Th' : [-0.25, 1.5],
      'K': [3.6,5],
      'Fracture porosity':[-8,-5],
      'Pyrite':[-6,-2],
      'Water-rock mass ratio':[0.0,0.004]
    }

    scaledict = {
      'neighbor porosity' : 'lin',
      'bulk density' : 'lin',
      'U' : 'log',
      'Th' : 'log',
      'K': 'log',
      'Fracture porosity': 'log',
      'Pyrite':'log',
      'Water-rock mass ratio':'lin'
    }

    fullparams = {
      'U': 'Uranium [ppm]',
      'K': 'Potassium [ppm]',
      "neighbor porosity": 'Primary porosity (vol/vol)',
      'Th': 'Thorium [ppm]',
      'Water-rock mass ratio': 'Water-rock mass ratio (g/g)',
      'bulk density': 'Bulk density [g/ccm]',
      'Fracture porosity': 'Apparent secondary porosity (vol/vol)',
      'Pyrite': 'Pyrite concentration [%]'
    }


    def Boxplots(all_model, BH_model, params, filename,
      colors = ['#d7191c', '#fdae61', '#2b83ba'],
      xlimsdict = xlimsdict,
      scaledict = scaledict,
      fullparams = fullparams,
      skew=False
      ):


        fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8,10))
        axs = axs.flatten()

        fitter_kwargs = {"BH":"all", 'exclude_expts':['PP']}

        if type(all_model) == type(''):
            all_model = [all_model]

        for i, param in enumerate(params):

            fitter_kwargs['log_unit_adjust'] = 0
            if param == 'K':
                fitter_kwargs['log_unit_adjust'] = 4

            xlims = xlimsdict[param][0], xlimsdict[param][1]
            label = fullparams[param]

            # fist fit a line for all boreholes together
            f_allBH, d_allBH, xs_allBH, ys_allBH = histplotter.line(axs[i],
              param, all_model[0], xlims=xlims,
              scale=scaledict[param], fitter_kwargs=fitter_kwargs,
              line_kwargs={'c':'k', 'label':all_model[0]+' for all BH', 'lw':2})

            # if a second fit for all BH is requested, add that line too.
            if len(all_model) > 1:
                f_allBH, d_allBH, xs_allBH, ys_allBH = histplotter.line(axs[i],
                  param, all_model[1], xlims=xlims,
                  scale=scaledict[param], fitter_kwargs=fitter_kwargs,
                  line_kwargs={'c':'k', 'label':all_model[1]+' fits for all BH', 'ls':'dotted', 'lw':2})


            # iterate through boreholes for their individual fits
            for j, (BH, c) in enumerate(zip(['BH01', 'BH02', 'BH03'], colors)):

                thisBH_kwargs = {"BH":BH, 'exclude_expts':['PP'], 'log_unit_adjust': fitter_kwargs['log_unit_adjust']}

                _v = histplotter.line(axs[i], param, BH_model, xlims=xlims,
                  scale=scaledict[param], fitter_kwargs=thisBH_kwargs,
                  line_kwargs={'c':c, 'label': BH_model+' for '+BH, 'lw':2})

            # add boxplots for the cumulative BH data.
            a2 = axs[i].twinx()
            histplotter.boxplots(a2, param, BH_model, xlims=xlims,
              scale=scaledict[param], colors=colors,
              fitter_kwargs=fitter_kwargs, hist_kwargs={'bins':10, 'alpha':0.3}, skew=skew)

            a2.set_ylim(-5,5)
            axs[i].set_title(label)
            axs[i].set_xlabel(label)
            axs[i].set_xlim(xlims)
            if scaledict[param] == 'log':
                axs[i].set_xlabel('log$_{10}$('+label+r')')


        for i, ax in enumerate(axs):
            ax.set_ylabel('Probability density')
            ax.text(0.01,0.99, string.ascii_lowercase[i], va='top', ha='left', fontsize=12, fontweight='bold', transform=ax.transAxes)
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*2)


        # The only specificity of the code is when plotting the legend
        h, l = np.hstack([axis.get_legend_handles_labels()
                          for axis in axs[0].figure.axes
                          if axis.bbox.bounds == axs[0].bbox.bounds]).tolist()
        print(h, l)
        l.insert(2, l[-1])
        l.pop(-1)
        h.insert(2, h[-1])
        h.pop(-1)
        axs[0].legend(handles=h, labels=l, bbox_to_anchor=[1.1, 1.3], loc='center', ncols=3)
        plt.tight_layout()
        fig.subplots_adjust(wspace=0.23, hspace=0.33)

        plt.savefig(os.path.dirname(__file__)+'/../../Figures/'+filename)
        plt.close()




    def input_data_basic_properties(params,
      scaledict = scaledict,
      fullparams = fullparams):

        BHlist = ['BH01', 'BH02', 'BH03', 'all']

        allprops = {}
        for BH in BHlist:
            hf = histfitter(BH = BH, exclude_expts=['PP'])
            for i, p in enumerate(params):
                allprops[fullparams[p]+' '+BH] = hf.characterize(p, scaledict[p])

        _df = pd.DataFrame.from_dict(allprops, orient='index')
        # _df = _df.sort_values('BH', ascending=True)
        _df = _df.sort_values(['param', 'BH'], ascending=(True, True))
        _df.to_csv(os.path.dirname(__file__)+'/../../HardRockMC/data/DataInput/basic_stats_properties.csv',
          columns = ['param', 'BH', 'N', 'mean', 'median', 'min', 'pc5', 'pc25', 'pc75', 'pc95', 'max'])




    def KS_test(params,
      xlimsdict = xlimsdict,
      scaledict = scaledict,
      fullparams = fullparams
      ):
        """
        Perform KS tests on normal distribution fits to the input data listed in
        the params argument. Saves KS p values to a csv in HardRockMC/data/DataInput.
        """

        all_model = ['norm']
        BH_model = 'norm'

        # it is convenient to use histplotter to get the data we want,
        # so make a plot anyway (it will later be closed)
        fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(8,10))
        axs = axs.flatten()

        fitter_kwargs = {"BH":"all", 'exclude_expts':['PP']}

        ksp_all = []
        ksp_BHall = [[],[],[]]
        ksp_BHBH = [[],[],[]]


        for i, param in enumerate(params):

            # line for all boreholes together
            f_allBH, d_allBH, xs_allBH, ys_allBH = histplotter.line(axs[i],
              param, all_model[0],
              scale=scaledict[param], fitter_kwargs=fitter_kwargs,
              line_kwargs={'c':'k', 'label':all_model[0]+' for all BH', 'lw':2})

            dist = eval("spys.norm")
            n = dist(f_allBH[0], f_allBH[1])
            ks = spys.kstest(d_allBH, n.cdf)
            print(param, ks.statistic, ks.pvalue)
            ksp_all.append(ks.pvalue)

            for j, BH in enumerate(['BH01', 'BH02', 'BH03']):

                thisBH_kwargs = {"BH":BH, 'exclude_expts':['PP']}

                _v = histplotter.line(axs[i], param, BH_model,
                  scale=scaledict[param], fitter_kwargs=thisBH_kwargs,
                  line_kwargs={'label': BH_model+' for '+BH, 'lw':2})

                ks = spys.kstest(_v[1], n.cdf)
                ksp_BHall[j].append(ks.pvalue)
                _n = dist(_v[0][0], _v[0][1])
                ks = spys.kstest(_v[1], _n.cdf)
                ksp_BHBH[j].append(ks.pvalue)

        plt.close()
        # convert to dict with labels etc
        ks_df = pd.DataFrame({'Parameter': params,
        'all data : norm all data': ksp_all,
        'BH01 : norm all data': ksp_BHall[0],
        'BH02 : norm all data' : ksp_BHall[1],
        'BH03 : norm all data' : ksp_BHall[2],
        'BH01 : norm BH01': ksp_BHBH[0],
        'BH02 : norm BH02' : ksp_BHBH[1],
        'BH03 : norm BH03' : ksp_BHBH[2],
        })

        ks_df.to_csv(os.path.dirname(__file__)+'/../../HardRockMC/data/DataInput/KS_p_values.csv')
