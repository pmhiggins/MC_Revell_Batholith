import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../')
sys.path.append(os.path.dirname(__file__)+'/../../')

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as spys
import pandas as pd

from histfitter import histfitter
from histplotter import histplotter
from RetrieveProductionRates import RetrieveProductionRates

from uncertainties.unumpy import uarray


class gas_flux:


    default_params = {
      'Y_H2': r'H$_2$ production rate [mol (m$^{3}$ rock yr)$^{-1}$]',
      'Y_4He': r'$^4$He production rate [mol (m$^{3}$ rock yr)$^{-1}$]',
      "Y_40Ar": r'$^{40}$He production rate [mol (m$^{3}$ rock yr)$^{-1}$]',
      'Y_SO4':  r'Sulfate production rate [mol (m$^{3}$ rock yr)$^{-1}$]',
      'Y_H2:Y_He': r'Ratio between productino rate of H$_2$ and He',
    }

    default_xlims = {'Y_H2': [-10,-7.8],
    'Y_4He': [-11.5,-8],
    "Y_40Ar": [-12,-10],
    'Y_SO4': [-16,-8.5],
    'Y_H2:Y_He':[0,5]
    }

    def get_CIS_lower_upper(fn, param_label, xlims, xdim):
        """
        Get two arrays of length CIS_num which represent the upper and lower
        95% confidence limit in the gas production probability density functions.
        """

        df = RetrieveProductionRates.ProductionRate_df('CIS/'+fn)

        lower_upper = {}

        df[param_label] = uarray(df[param_label+' nom'], df[param_label+' sigma'])

        xs = np.linspace(xlims[param_label][0], xlims[param_label][1], xdim)

        yarr = np.zeros((len(df[param_label]), len(xs)))

        for i, ms in enumerate(df[param_label]):

            yarr[i] = spys.norm.pdf(xs, ms.n, ms.s)

        yarr = np.transpose(yarr)
        smallest = []
        largest = []
        for x, y in zip(xs, yarr):
            smallest.append(np.percentile(y, 2.5))
            largest.append(np.percentile(y,  97.5))

        return (smallest, largest)


    def individual_output_boxes(fn_KDE='', fn_norm='', fn_CIS='', weighted=False,
      params = default_params,
      xlims = default_xlims,
      CIS_xdim=1000,
      fig_preamble='',
      scale='log'):

        """
        Plot and save in independent figures the outputs, using data saved in
        fn_KDE for the KDE distributed inputs and fn_norm for the norm
        distributed inputs.

        Generates probability densities and boxplots showing key quantile
        ranges in a single-axis matplotlib figure.

        Optionally, confidence interval based on the normal fits to inputs
        can be plotted as a fill_between bu passing fn_CIS.
        """

        save_dirn = os.path.dirname(__file__)+'/../../Figures/'

        for param in params.keys():

            fig, axs = plt.subplots(nrows=1,ncols=1, figsize=(6,6))

            if weighted:

                # first plot the lines form the main MC results
                histplotter.line(axs, param, 'KDE', scale=scale,
                  xlims = xlims[param],
                  fitter_kwargs={'outputs': 'nonweighted/'+fn_KDE, 'BH':'all'},
                  line_kwargs={'label':'KDE inputs; \n'+'all BH pooled', 'c':'k'})

                histplotter.line(axs, param, 'norm', scale=scale,
                  xlims = xlims[param],
                  fitter_kwargs={'outputs': 'nonweighted/'+fn_norm, 'BH':'all'},
                  line_kwargs={'label':'Nominal normal inputs;\n'+'all BH pooled', 'c':'b'})

                # now the 'weighted' ones
                histplotter.line(axs, param, 'KDE', scale=scale,
                  xlims = xlims[param],
                  fitter_kwargs={'outputs': 'weighted/'+fn_KDE, 'BH':'weighted'},
                  line_kwargs={'label':'KDE inputs; \n'+'equal BH weighting', 'c':'k', 'ls':'dotted'})

                histplotter.line(axs, param, 'norm', scale=scale,
                  xlims = xlims[param],
                  fitter_kwargs={'outputs': 'weighted/'+fn_norm, 'BH':'all'},
                  line_kwargs={'label':'Nominal normal inputs;\n'+'equal BH weighting', 'c':'b', 'ls':'dotted'})

            else:

                # only plot the lines form the main MC results
                histplotter.line(axs, param, 'KDE', scale=scale,
                  xlims = xlims[param],
                  fitter_kwargs={'outputs': 'nonweighted/'+fn_KDE, 'BH':'all'},
                  line_kwargs={'label':'KDE inputs', 'c':'k', 'lw':'2'})

                histplotter.line(axs, param, 'norm', scale=scale,
                  xlims = xlims[param],
                  fitter_kwargs={'outputs': 'nonweighted/'+fn_norm, 'BH':'all'},
                  line_kwargs={'label':'Normal inputs', 'c':'b', 'lw':'2'})


            # plot a confidence interval if needed.
            # this will noticably increase compute time.
            if fn_CIS != '':

                lower_upper = gas_flux.get_CIS_lower_upper(fn_CIS, param, xlims, CIS_xdim)
                # now we need to find the limits for the fillbetween norm confidence intervals
                xs = np.linspace(xlims[param][0], xlims[param][1], CIS_xdim)

                axs.fill_between(xs, lower_upper[0], lower_upper[1], facecolor='tab:blue', alpha=0.5,
                  label='95% CI in \n'+'normal approx. \n'+'due to sample size')

            if scale=='log':
                axs.set_xlabel(r'$\log_{10}$('+params[param]+')', fontsize=14)
            else:
                axs.set_xlabel(params[param])

            # add boxplots for the cumulative BH data.
            a2 = axs.twinx()
            histplotter.boxplots(a2, param, ['KDE', 'norm'], BHs=['all', 'all'], xlims=xlims[param],
              scale='log', colors=['k', 'b'],
               fitter_kwargs=[{'outputs': 'nonweighted/'+fn_KDE, 'BH':'all'}, {'outputs': 'nonweighted/'+fn_norm, 'BH':'all'}], extended=True)

            a2.set_ylim(-6,4)

            # rename the legend from the histplotter defaults
            h, l = np.hstack([axis.get_legend_handles_labels()
                              for axis in axs.figure.axes
                              if axis.bbox.bounds == axs.bbox.bounds]).tolist()

            l[3] = 'KDE result boxplot'
            l[4] = 'norm result boxplot'

            axs.set_xlim(xlims[param][0], xlims[param][1])
            axs.set_ylim(axs.get_ylim()[0], axs.get_ylim()[1]*1.5)
            axs.set_ylabel('Probability Density', fontsize=14)
            axs.set_title(params[param], fontsize=14)
            legend =axs.legend(handles=h, labels=l, bbox_to_anchor=[0.5, 1.2], loc='center', ncols=2)

            plt.tight_layout()
            fig.subplots_adjust(top=0.75, bottom=0.1)
            plt.savefig(save_dirn+'/'+fig_preamble+param+'_'+scale+'.pdf')
            # plt.savefig(save_dirn+'/'+param+'_'+scale+'.png', dpi=300)
            plt.close()


    def output_data_basic_properties(params, fn_KDE, fn_norm,
      full_params=default_params):
        """
        Prepare a DataFrame of basic statistical properties in the lists params.
        Save the output and print it to terminal.
        """

        allprops = {}
        for fn, fit in zip([fn_KDE, fn_norm], ['KDE', 'norm']):
            for weight in ['nonweighted','weighted']:
                hf = histfitter(BH = 'all', exclude_expts=['PP'], outputs = weight+'/'+fn)
                for i, p in enumerate(params):
                    allprops[full_params[p]+' '+fit+weight] = hf.characterize(p, 'log')

        _df = pd.DataFrame.from_dict(allprops, orient='index')
        _df = _df.sort_values('param', ascending=True)
        _df.to_csv(os.path.dirname(__file__)+'/../../HardRockMC/data/DataOutput/basic_stats_properties.csv',
          columns = ['param', 'N', 'mean', 'median', 'min', 'pc5', 'pc25', 'pc75', 'pc95', 'max'])
        print(_df)
