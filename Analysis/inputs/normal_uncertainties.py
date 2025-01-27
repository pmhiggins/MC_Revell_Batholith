import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.dirname(__file__)+'/../../')

import numpy as np
import scipy.stats as spys
import matplotlib.pyplot as plt
import math

from histfitter import histfitter



class normal_uncertainties:

    xlimsdict = {
      'neighbor porosity' : [0.0, 0.01],
      'bulk density' : [2.55,2.75],
      'U' : [-0.5, 1.5],
      'Th' : [-0.25, 1.5],
      'K': [-0.5,math.log10(20)],
      'Fracture porosity':[-8,-5],
      'Pyrite': [-7,-2]
    }

    scaledict = {
      'neighbor porosity' : 'lin',
      'bulk density' : 'lin',
      'U' : 'log',
      'Th' : 'log',
      'K': 'log',
      'Fracture porosity': 'log',
      'Pyrite':'log'
    }

    fullparams = {
      'U': r'log$_{10}$ Uranium [ppm]',
      'K': r'log$_{10}$ Potassium [ppm]',
      "neighbor porosity": 'Primary porosity',
      'Th': r'log$_{10}$ Thorium [ppm]',
      'Water-rock mass ratio': 'Water-rock mass ratio',
      'bulk density': 'Bulk density [g/ccm]',
      'Fracture porosity': 'Apparent bulk porosity',
      'Pyrite':r'log$_{10}$ Concentration of Pyrite'
    }


    params = ['neighbor porosity',
      'bulk density',
      'U',
      'Th',
      'K']

    def normal_confidencelimits(
      small_n = 10,
      large_n = 500,
      params=params,
      xlimsdict=xlimsdict,
      scaledict = scaledict,
      fullparams=fullparams,
      CIS_num=100000,
      CIS_xdim=1000
    ):

        for param in params:

            # lists for each line: BH1, BH2, BH3, total
            df_small = [small_n,small_n,small_n,small_n*3]
            df_large = [large_n,large_n,large_n,large_n*3]
            c = ['#d7191c', '#fdae61', '#2b83ba', 'k']
            BH = ['BH01', 'BH02', 'BH03', 'all']
            # label = ['BH01; n=80', 'BH02; n=25', 'BH03; n=20', 'total; n=125']

            # xs = np.linspace(xlimsdict[param][0], xlimsdict[param][1], num=int(CIS_xdim))

            fig, axs = plt.subplots(ncols=3, figsize=(8,3), sharey=True)

            # i refers to BH1, BH2, BH3, total
            for i in range(4):

                hf = histfitter(BH[i], only_expt=None, exclude_expts=['PP'],
                  outputs = None)

                density, param_values, xs, _pdf = hf.fit(
                  param, 'norm', xlims = xlimsdict[param],
                  scale=scaledict[param],
                  min_infs_converge='remove', xdim=CIS_xdim)

                df_true = len(param_values)

                mean = density[0]
                std = density[1]

                for df, ax in zip([df_small[i], df_true, df_large[i]], axs):

                    nom_pdf = spys.norm.pdf(xs, mean, std)
                    ax.errorbar([mean], [np.max(nom_pdf)], xerr=[std/math.sqrt(df)], ecolor=c[i], capsize=5)
                    ax.plot(xs, spys.norm.pdf(xs, mean, std), c=c[i])

                    # build distributions to assert
                    # uncertainties in mean and standard deviation
                    s_dist = np.sqrt((df-1)/spys.chi2.rvs(df, size=int(CIS_num)))*std
                    m_dist = spys.norm.rvs(mean, (std/(math.sqrt(df))), size=CIS_num)

                    # draw a distribution using each selection of m and s
                    yarr = np.zeros((CIS_num, CIS_xdim))
                    for j, (m,s) in enumerate(zip(m_dist, s_dist)):
                        yarr[j] = spys.norm.pdf(xs, m, s)


                    # take transpose to find the max/min at given x
                    yarr = np.transpose(yarr)

                    smallest = []
                    largest = []

                    for x, y in zip(xs, yarr):
                        smallest.append(np.percentile(y, 2.5))
                        largest.append(np.percentile(y,  97.5))

                    lower_upper = [smallest, largest]

                    ax.fill_between(xs, smallest, largest, facecolor=c[i], alpha=0.3,
                      label=r'n ='+str(df))

                    ax.set_xlabel(fullparams[param])
                    ax.set_ylabel('Probability density')
                    if ax == axs[1]:
                        ax.legend(fontsize='smaller')

            axs[0].set_title('Dataset size: '+str(df_small[0])+' per BH')
            axs[1].set_title('Actual raw dataset size')
            axs[2].set_title('Dataset size: '+str(df_large[0])+' per BH')


            plt.tight_layout()
            plt.subplots_adjust(wspace=0.2, bottom=0.165)
            savedir = os.path.dirname(__file__)+'/../../Figures/'
            plt.savefig(savedir+'S5_input_sampleconfidence_'+param+'.pdf')
            plt.close()
            print(param)
