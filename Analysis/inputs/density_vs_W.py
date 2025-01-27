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

class density_vs_W:


    def dens_vs_W():

        KThU_df, D_df, W_df = EPM.get_dataframes()

        _d = 'neighbor'

        # fig for all three BH, fig2 for individual boreholes
        fig, axs = plt.subplots(ncols=2, figsize=(11,5))
        fig2, axs2 = plt.subplots(ncols=3, figsize=(12,4))

        markers = ['X', 's', '^', 'o']
        colors = ['blue', 'tab:orange', 'green']
        for j, BH in enumerate(['BH01', 'BH02', 'BH03']):

            for i, exp in enumerate(['IsoEx LAB', 'PWoutdiff', 'PW', 'AQ']):
                _W_df = W_df[(W_df['Borehole'] == BH) & ((W_df['Type'] == exp))]

                # all boreholes in one plot
                axs[0].errorbar(_W_df['Water-rock mass ratio nom'], _W_df[_d+' bulk density nom'], yerr=_W_df[_d+' bulk density sigma'], xerr=_W_df['Water-rock mass ratio sigma'], fmt=markers[i], capsize=5, c=colors[j])
                axs[0].set_title('With errorbars')

                # single borehole per plot
                axs2[j].errorbar(_W_df['Water-rock mass ratio nom'], _W_df[_d+' bulk density nom'], yerr=_W_df[_d+' bulk density sigma'], xerr=_W_df['Water-rock mass ratio sigma'], fmt=markers[i], capsize=5, c=colors[j])
                axs2[j].set_title(BH)

                if j != 0:
                    axs[1].scatter(_W_df['Water-rock mass ratio nom'], _W_df[_d+' bulk density nom'],
                      marker=markers[i], c=colors[j])
                else:
                    axs[1].scatter(_W_df['Water-rock mass ratio nom'], _W_df[_d+' bulk density nom'],
                      marker=markers[i], label=exp, c=colors[j])
                axs[1].set_title('Without errorbars')
            axs[1].plot([np.nan], [np.nan], marker='1', c=colors[j], label=BH)

        axs[1].legend()
        for ax in chain(axs, axs2):
            ax.set_xlabel(r'Water content [$m_w$ / $m_r$]')
            ax.set_ylabel('Bulk density [g/ccm]')
            ax.set_xlim(axs[0].get_xlim())
            ax.set_ylim(axs[0].get_ylim())

        fig.suptitle('Bulk density versus water content')

        fig.tight_layout()
        fig2.tight_layout()

        fig.savefig(os.path.dirname(__file__)+'/../../Figures/S4a_dens_vs_W.pdf')
        fig2.savefig(os.path.dirname(__file__)+'/../../Figures/S4b_dens_vs_W.pdf')
