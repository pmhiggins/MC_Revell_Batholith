
import os
import ExceltoPandasMethods as EPM
from histfitter import histfitter
from histplotter import histplotter
from copy import deepcopy
from uncertainties import ufloat as uf
import scipy.stats as spys
import math
import numpy as np



class SampleGenerator:


    default_params = [ 'neighbor porosity',
     'bulk density',
     'U',
     'Th',
     'K',
     'Water-rock mass ratio',
     'Pyrite'
    ]


    all_lin = ['lin',]*len(default_params)


    def __init__(self,
      num, fitter,
      exclude_expts=['PP'],
      BH = 'all',
      params = default_params,
      scales = all_lin):

        self.num = num
        self.fitter = fitter
        self.exclude_expts = exclude_expts
        self.BH = BH
        self.params = params
        self.scales = scales

        # if fitter is a string, all are to be fit with the same distribution
        # so make a list reflecting that
        if type(self.fitter) == type(''):
            self.fitter = [self.fitter,]*len(self.params)


    def Generate_nonweighted(self, fn):
        """
        Generate an input sample for the gas production by pooling all borehole
        data together and NOT weighting for equal contributions for each borehole.
        """

        samples = []

        saveas = os.path.dirname('HardRockMC/data/DataInput/Redistributed/')+'/'+fn

        for _param, _fitter, _scale in zip(self.params, self.fitter, self.scales):

            _sample = []

            hf = histfitter(self.BH, exclude_expts=self.exclude_expts)

            _sample = hf.resample(self.num, _param, _fitter,
              scale=_scale,
              fitter_kwargs={})

            if _scale == 'log':
                samples.append(10**_sample)
            else:
                samples.append(_sample)

        self.params.append('PyriteSA')
        samples.append(np.random.uniform(low=5.65, high=28.3, size=len(_sample)))

        print(saveas)
        _df = EPM.distributions_to_csv(self.params, samples,
          filename=saveas, negpurge=True)

        return _df

    def Generate_weighted(self, fn):
        """
        Generate an input sample for the gas production by pooling all borehole
        data together and weighting for equal contributions for each borehole.
        """

        samples = []

        saveas = os.path.dirname('HardRockMC/data/DataInput/Redistributed_weighted/')+'/'+fn

        for _param, _fitter, _scale in zip(self.params, self.fitter, self.scales):

            _sample = []

            for _BH in ['BH01', 'BH02', 'BH03']:

                hf = histfitter(_BH, exclude_expts=self.exclude_expts)
                __sample = hf.resample(int(self.num/3), _param, _fitter,
                  scale=_scale,
                  fitter_kwargs={})
                _sample.extend(__sample)

            _sample = np.array(_sample)

            if _scale == 'log':
                samples.append(10**_sample)
            else:
                samples.append(_sample)

        self.params.append('PyriteSA')
        samples.append(np.random.uniform(low=5.65, high=28.3, size=len(_sample)))

        _df = EPM.distributions_to_csv(self.params, samples,
          filename=saveas, negpurge=True)

        return _df


    def Generate_CIS(self, fn):
        """
        Generate an input sample for the gas production by pooling all borehole
        data together. This generates a distribtution in mean values and a
        distribution in standrad deviation values as described in the
        Supplemental Material in Higgins, Song et al., (2025).
        """

        samples = []
        new_params = []

        saveas = os.path.dirname('HardRockMC/data/DataInput/Redistributed_CIS/')+'/'+fn

        for _param, _fitter, _scale in zip(self.params, self.fitter, self.scales):

            this_sample = []
            hf = histfitter(self.BH, exclude_expts=self.exclude_expts)

            density, param_values, xs, _pdf = hf.fit(
              _param, _fitter,
              scale=_scale,
              min_infs_converge='remove')

            mean = density[0]
            std = density[1]
            df = len(param_values) # degrees of freedom

            # distribution of possible standard deviations
            s_dist = np.sqrt((df-1)/spys.chi2.rvs(df, size=int(self.num)))*std
            # distribution of possible mean values
            m_dist = spys.norm.rvs(mean, (std/(math.sqrt(df))), size=self.num)
            new_params.append(_param+ ' nom')
            new_params.append(_param+ ' sigma')
            samples.append(m_dist)
            samples.append(s_dist)

        _df = EPM.distributions_to_csv(new_params, samples,
          filename=saveas, negpurge=False)

        return _df


    def GenerateWarr2023Samples(self, fn, type='norm'):
        """
        Generate samples for input parameters using the statistical moments
        for Kidd Creek (Gaussians) described in Warr et al., (2023).
        """

        samples = []

        fn = 'Warr2023_'+type+'_'+str(self.num)+'.csv'

        saveas = os.path.dirname('data/DataInput/Redistributed')+'/'+fn

        for param in self.params:
            print(param)
            _mu = histplotter.Warr23_means_stds[param][0]
            _sigma = histplotter.Warr23_means_stds[param][1]

            if type=='norm':
                samples.append(np.random.normal(_mu, _sigma, self.num))
            elif type=='uniform':
                _min = _mu - (3*_sigma)
                _max = _mu + (3*_sigma)
                if param == 'Pyrite' or param =='PyriteSA':
                    # for these two warr2023 used a 5 sigma approach
                    # where the upper and lower limits were assumed to be
                    # at the 5 sigma level.
                    _min = _mu - (5*_sigma)
                    _max = _mu + (5*_sigma)
                samples.append(np.random.uniform(_min, _max, self.num))
            else:
                raise ValueError('unknown distribution passed while attempting to recreate Warr et al 2023 results.')

        self.params.append('Water-rock mass ratio')
        self.params.append('rock density')

        bulk_dens = samples[1]*(1-samples[0]) + (samples[0]*samples[5])
        W = samples[0]*samples[5] / ((1-samples[0])*samples[1])
        samples.append(W)
        samples.append(deepcopy(samples[1])) #rock density
        samples[1] = bulk_dens

        EPM.distributions_to_csv(self.params, samples,
          filename=saveas, negpurge=True)

        return fn
