
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mplc
import scipy
from copy import deepcopy
from uncertainties import ufloat as uf
from uncertainties import unumpy as unp

from scipy import stats as spys
import seaborn as sns
from itertools import chain
import numpy as np
import math

from histfitter import histfitter
import ExceltoPandasMethods as EPM


class histplotter:

    Warr23_means_stds = {
      'neighbor porosity' : [0.01, 0.0045/3.],
      'bulk density' : [3.0,0.3/3.],
      'U' : [1.455, 0.545/3.],
      'Th' : [6.655, 2.345/3.],
      'K': [17.4e3,(2.600e3)/3.], # converted to ppm
      'water density': [1.11, 0.11/3.],
      'Pyrite':[20.75, 17.25/5.],
      'PyriteSA':[23.6915, 23.4085/5.],
      'Y_SO4':[5e-10, 2e-10],
      'Y_4He':[4.9e-11,0.4e-11],
      'Y_H2':[3.5e-9,0.6e-9],
      'Y_40Ar':[8.8e-12, 0.5e-12]
    }


    def line(ax, param, model, scale='lin', xlims=None,
      fitter_kwargs={}, line_kwargs={}):

        hf = histfitter(**fitter_kwargs)
        _, v, xs, ys = hf.fit(param, model, scale=scale, xlims=xlims, min_infs_converge='remove')

        ax.plot(xs, ys, **line_kwargs)

        return _, v, xs, ys

    def hist(ax, param, scale='lin', xlims=None, bins=10,
      fitter_kwargs={}, hist_kwargs={}):
        """
        Plot a histogram of param. Borehole (or all) should be passed as a
        fitter_kwarg.

        If on a log scale and -infs are encounterd (<=0), add the count
        of them as text in the upper right hand corner.
        """
        hf = histfitter(**fitter_kwargs)
        v, num_infs = hf.get_param_values(param, scale=scale, min_infs_converge='remove', return_num_infs=True)

        # _, v, xs, ys = hf.fit(param, model, xlims=xlims, scale=scale)
        ax.hist(v, density=True, **hist_kwargs)
        if num_infs >0 :
            pc_infs = 100* num_infs / (len(v)+num_infs)

            ax.text(1,1, str(round(pc_infs,3))+'% less than zero',
              ha='right', va='top', transform=ax.transAxes)

        return v


    def stacked_hist(ax, param, model, BHs=['BH01', 'BH02', 'BH03'],
      scale='lin', xlims=None, bins=10, colors= ['blue', 'tab:orange', 'green'],
      fitter_kwargs={}, hist_kwargs={}):


        histlist = []
        _fitter_kwargs = deepcopy(fitter_kwargs)

        for i, BH in enumerate(BHs):
            _fitter_kwargs['BH'] = BH
            hf = histfitter(**_fitter_kwargs)
            _, v, xs, ys = hf.fit(param, model, xlims=xlims, scale=scale)

            histlist.append(v)


        ax.hist(histlist, density=True, stacked=True, color=colors[:len(BHs)], label=BHs, **hist_kwargs)



    def violins(ax, param, model, Warr=True, BHs=['BH01', 'BH02', 'BH03'],
      scale='lin', xlims=None, bins=10, colors= ['blue', 'tab:orange', 'green'],
      fitter_kwargs={}, hist_kwargs={}):

        histlist = []
        _fitter_kwargs = deepcopy(fitter_kwargs)

        Warr_norm = []
        try:
            Warr_norm = np.random.normal(
              histplotter.Warr23_means_stds[param][0],
              histplotter.Warr23_means_stds[param][1],
              10000
            )
            if scale=='lin':
                histlist.append(Warr_norm)
            elif scale == 'log':
                histlist.append(np.log10(Warr_norm))

        except:
            print('Missing Warr parameter: '+param)


        for i, BH in enumerate(BHs):
            _fitter_kwargs['BH'] = BH
            hf = histfitter(**_fitter_kwargs)
            _, v, xs, ys = hf.fit(param, model, xlims=xlims, scale=scale)

            histlist.append(v)

        parts = ax.violinplot(histlist, vert=False, showmeans=True, showextrema=True,
          # quantiles=[0.25,0.75])
        )

        if Warr_norm != []:
            parts['bodies'][0].set_facecolor('m')
            parts['bodies'][0].set_edgecolor('black')
            parts['bodies'][0].set_alpha(0.75)

            for pc, col in zip(parts['bodies'][1:], colors):
                pc.set_facecolor(col)
                pc.set_edgecolor('black')
                pc.set_alpha(0.75)

        else:
            for pc, col in zip(parts['bodies'], colors):
                pc.set_facecolor(col)
                pc.set_edgecolor('black')
                pc.set_alpha(0.75)


    def boxplots(ax, param, model, Warr=True, BHs=['BH01', 'BH02', 'BH03'],
      scale='lin', xlims=None, bins=10, colors= ['blue', 'tab:orange', 'green'], extended=False,
      fitter_kwargs={}, hist_kwargs={}, skew=False):

        histlist = []
        print(fitter_kwargs)
        Warr_norm = []
        labels = [BH + ' input data' for BH in BHs]

        for i, BH in enumerate(BHs):
            m=model
            _fitter_kwargs=deepcopy(fitter_kwargs)
            if type(model) == type([]):
                m=model[i]
                print(fitter_kwargs)
                _fitter_kwargs = fitter_kwargs[i]
            _fitter_kwargs['BH'] = BH
            hf = histfitter(**_fitter_kwargs)
            _, v, xs, ys = hf.fit(param, m, xlims=xlims, scale=scale)

            histlist.append(v)

        try:
            Warr_norm = np.random.normal(
              histplotter.Warr23_means_stds[param][0],
              histplotter.Warr23_means_stds[param][1],
              100000
            )
            if scale=='lin':
                histlist.append(Warr_norm)
            elif scale == 'log':
                histlist.append(np.log10(np.abs(Warr_norm)))
            colors.append('slategray')
            labels.append('Warr et al., (2023) \n Kidd Creek')
        except:
            print('Missing Warr parameter: '+param)


        symb='.'
        if extended:
            symb=''
            bg = ax.boxplot(histlist, sym=symb, vert=False, whis=(1,99), widths=0.75, meanline=True, showmeans=True, meanprops={'color':'w', 'linestyle':'dotted'}, whiskerprops={'linestyle':'dashed'})


        parts = ax.boxplot(histlist, sym=symb, vert=False, whis=(5,95), widths=0.75, labels=labels, meanline=True, showmeans=True, meanprops={'color':'w', 'linestyle':'dotted'})

        if not extended:
            parts['fliers'][-1].set(alpha=0)

        for i, (artist, med, color, label, hist_data) in enumerate(zip(parts['boxes'], parts['medians'], colors, labels, histlist)):
            patch = mpatches.PathPatch(artist.get_path(), color=mplc.to_rgba(color, 0.75),
            label=label)
            ax.add_artist(patch)
            med.set_color('k')
            if skew and i<3:
                S_KP = spys.skew(hist_data)
                print(param, 'BH: ', BH, S_KP)
                ax.text(xlims[1], i+1, str(round(S_KP,2))+'  ', ha='right', va='center', fontsize=10)
            elif skew:
                ax.text(xlims[1], i+1, r'$S_{KP}$   ', ha='right', va='center', fontsize=10)



        ax.set_yticks([])
        # for artist, color in zip(parts['boxes'], colors):
            # patch = mpatches.PathPatch(artist.get_path(), color=mplc.to_rgba(color, 0.75))
            # ax.add_artist(patch)



    def overlay_warr(ax, param, scale='lin', xlims=None, line_kwargs={}):
        mean, std = None, None
        try:
            mean = histplotter.Warr23_means_stds[param][0]
            std = histplotter.Warr23_means_stds[param][1]
        except:
            print('Missing Warr parameter: '+param)
            return

        if xlims == None:
            xlims = [mean-(6*std), mean+(6*std)]

        xs = np.linspace(xlims[0], xlims[1], num=1000)

        if scale == 'lin':
            ax.plot(xs, spys.norm.pdf(xs, mean, std), **line_kwargs)
        elif scale == 'log':
            ax.plot(np.log10(xs), spys.norm.pdf(xs, mean, std), **line_kwargs)


    def line_pc_difference(ax, params, models,
      scale='lin', xlims=None,
      fitter_kwargs1={}, fitter_kwargs2={}, line_kwargs={}):

        hf1 = histfitter(**fitter_kwargs1)
        hf2 = histfitter(**fitter_kwargs2)

        _1, v1, xs1, ys1 = hf1.fit(params[0], models[0], scale=scale, xlims=xlims)
        _2, v2, xs2, ys2 = hf2.fit(params[1], models[1], scale=scale, xlims=xlims)

        ax.plot(xs1, 100*(ys1 - ys2)/ys2, **line_kwargs)


    def scatter_samplescale(ax, xparam, yparam,
      errorbars=True,
      mean=False,
      fitter_kwargs={}, plot_kwargs={}):

        xparam_nom = xparam + ' nom'
        xparam_sigma = xparam + ' sigma'
        yparam_nom = yparam + ' nom'
        yparam_sigma = yparam + ' sigma'

        hf1 = histfitter(**fitter_kwargs)
        _df = hf1.which_df(xparam)

        try:
            if errorbars:
                ax.errorbar(_df[xparam_nom], _df[yparam_nom],
                  xerr=_df[xparam_sigma], yerr = _df[yparam_sigma], ls=plot_kwargs.pop('ls', 'none'),
                  capsize=plot_kwargs.pop('capsize', 3),
                  ms=plot_kwargs.pop('ms', 5), **plot_kwargs)
            else:
                ax.scatter(_df[xparam_nom], _df[yparam_nom],
                  s=plot_kwargs.pop('ms', 30), **plot_kwargs)

        except:
            ax.scatter(_df[xparam], _df[yparam],
              xerr=_df[xparam], yerr = _df[yparam],
              ms=plot_kwargs.pop('ms', 5), **plot_kwargs)

        if mean:
            return _df[xparam_nom].mean(), _df[xparam_nom].std(), _df[yparam_nom].mean(), _df[yparam_nom].std()
        else:
            return _df[xparam_nom], None, _df[yparam_nom], None



    def regression_samplescale(ax, xparam, yparam,
      fitter_kwargs={}, plot_kwargs={}):

        xparam_nom = xparam + ' nom'
        yparam_nom = yparam + ' nom'

        hf1 = histfitter(**fitter_kwargs)
        _df = hf1.which_df(xparam)

        g = sns.regplot(x=xparam_nom, y=yparam_nom, data=_df, marker='none', ci=95, ax=ax, **plot_kwargs)

        return g
