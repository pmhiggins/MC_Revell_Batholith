import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../../')


import HardRockMC.ExceltoPandasMethods as EPM
import matplotlib.pyplot as plt

from scipy import stats as spys
import seaborn as sns
from itertools import chain
import numpy as np
import math

from matplotlib import rcParams


class depth_scatter:


    def get_means(_df, yparam, xparam, interval):

        y_means = []
        y_stds = []
        x_means =[]
        x_stds =[]


        xmax = _df[xparam].max()*1.0001
        xmin = _df[xparam].min()*0.9999
        fullmean = _df[yparam].mean()
        fullmedian = _df[yparam].median()


        this_bottomx = xmin
        this_topx = xmin
        while this_topx < xmax:
            this_topx += interval
            __df = _df[(_df[xparam] >= this_bottomx) & (_df[xparam] <= this_topx)]
            y_means.append(__df[yparam].mean())
            y_stds.append(__df[yparam].std())
            x_means.append(__df[xparam].mean())
            x_stds.append(__df[xparam].std())
            this_bottomx += interval

        return y_means, y_stds, x_means, x_stds, fullmean, fullmedian




    def plot_depth_ind_joints(means=False, intervals=100):
        """ seaborn jointplots and pearson regression for individual input parameters"""

        KThU_df, D_df, W_df = EPM.get_dataframes(exclude_expts=['PP'])

        params = ['K', 'neighbor bulk density', 'Th', 'Water-rock mass ratio', 'U', 'neighbor porosity', 'Pyrite']
        labels = ['K [%]', 'Bulk Density [g/cm$^{3}$]', 'Th [ppm]', 'W [m$_w$ / m$_r$]', 'U [ppm]', 'Primary porosity [V$_w$ / V$_{bulk}$]', 'Pyrite [%]']

        df = [KThU_df, W_df, KThU_df, W_df, KThU_df, W_df, KThU_df]

        for i, key in enumerate(params):
            y = None
            try:
                y = key+' nom'
                g = sns.jointplot(data = df[i], x='true depth', y=y, kind="reg")
            except:
                y = key
                g = sns.jointplot(data = df[i], x='true depth', y=y, kind="reg")
            Pearson_r = spys.pearsonr(df[i]['true depth'], df[i][y])
            g.ax_joint.text(1,1,
              'Pearson: r=' +str(Pearson_r.statistic.round(3))+'; p='+str(Pearson_r.pvalue.round(3)),
              va='top', ha='right', transform=g.ax_joint.transAxes)

            if means:
                means,stds, depth_means, depth_stds, mean, med = depth_scatter.get_means(df[i], y, 'true depth', intervals)
                g.ax_joint.errorbar(depth_means, means, xerr=depth_stds, yerr=stds, c='r', capsize=5, fmt='x')

                g.ax_joint.axhline(mean, c='r')
                g.ax_joint.axhline(med, c='r', ls='dashed')

            g.fig.set_figwidth(4)
            g.fig.set_figheight(4)
            plt.xlabel('vertical depth [m]')
            plt.ylabel(labels[i])
            plt.tight_layout()
            plt.savefig(os.path.dirname(__file__)+'/../../Figures/depth_scatter/'+params[i].replace(' ', '')+'.pdf')
