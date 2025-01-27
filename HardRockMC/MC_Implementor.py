import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import ExceltoPandasMethods as EPM
from SampleGenerator import SampleGenerator
from GasProduction import GasProduction
import numpy as np
import scipy.stats as spys
import pandas as pd

class MC_Implementor:

    # name of parameters as in the source data file
    default_params = [ 'neighbor porosity',
     'bulk density',
     'U',
     'Th',
     'K',
     'Water-rock mass ratio',
     'Pyrite'
    ]

    # how to scale the params (linear or logarithmic) for drawing distributions
    # must have same shape as default_params
    default_scales = ['lin', 'lin', 'log','log','log','lin','log']

    def __init__(self,
      exclude_expts=['PP'],
      BH = 'all',
      params = default_params,
      scales = default_scales):

        self.params = params
        self.scales = scales
        self.exclude_expts = exclude_expts
        self.BH = BH


    def implement(self, num, fitter, fn, SampleType, porosity='primary'):
        """
        Perform and save an MC simulation for gas production using samples
        of size `num` and the `fitter` distribution as a best fit for the
        input parameters.

        SampleType must be one of three options and determines what kind of
        sample is taken and the dataframe that results. It must be one of
        "weighted", "nonweighted", "CIS".
        """

        if SampleType == 'CIS':
            return self.implement_CIS(num, fitter, fn)

        input_df = self.get_Sample(num, fitter, fn, SampleType)

        if porosity == 'secondary':
            input_df['Water-rock mass ratio'] = GasProduction.calc_W(input_df['bulk density'], input_df['Fracture porosity'])

        # Create a GasProduction object
        GP = GasProduction(input_df['bulk density'], input_df['Water-rock mass ratio'], input_df['U'], input_df['Th'], input_df['K'], input_df['PyriteSA'], input_df['Pyrite'],
          density_unit='ccm',
          K_unit = '%',
          sulfide_massfrac_unit='g/g',
          Energies=None, StoppingPower=None, Gs=None)

        # this calucates the gas production rates and adds them to its internal DataFrame
        GP.output_df()

        # save the output to a new csv of the same name,
        # but now in directory data/DataProcessed
        GP.save_output_df(fn, preamble='HardRockMC/data/DataOutput/'+SampleType+'/')

        return GP.MCdf



    def get_Sample(self, num, fitter, fn, SampleType):
        """
        Create a sample with size `num` of model input values. Rawdata is read
        in, a fit of type `fitter` is drawn, then a new sample is drawn from
        that fitted distribution.

        SampleType must be one of three options and determines what kind of
        sample is taken and the dataframe that results. It must be one of
        "weighted", "nonweighted", "CIS".
        """

        SG = SampleGenerator(
              num, fitter,
              exclude_expts=self.exclude_expts,
              BH = self.BH,
              params = self.params,
              scales = self.scales)

        if SampleType == 'nonweighted':
            return SG.Generate_nonweighted(fn)
        elif SampleType == 'weighted':
            return SG.Generate_weighted(fn)
        elif SampleType == 'CIS':
            return SG.Generate_CIS(fn)
        else:
            raise ValueError('Unknown input sample type requested!')



    def implement_CIS(self, num, fitter, fn, internal_num=1000000):

        input_df = self.get_Sample(internal_num, 'norm', fn, 'CIS')
        # this is a dataframe with num rows, and each entry is a
        # list of ufloats, in linear space or log space, depending on the
        # self.scales parameters.

        out_gases = ['Y_H2', 'Y_4He', 'Y_SO4','Y_40Ar']
        out_gases_dict = {'Y_H2 nom':[], 'Y_H2 sigma':[],
          'Y_4He nom':[],'Y_4He sigma':[],
          'Y_SO4 nom':[], 'Y_SO4 sigma':[],
          'Y_40Ar nom':[], 'Y_40Ar sigma':[]}

        for i in range(internal_num):
            # now we want to do an MC in each sample entry
            _internalMC = {}
            for param, scale in zip(self.params, self.scales):
                _mu = input_df[param+' nom'][i]
                _sigma = input_df[param+' sigma'][i]

                _internalMC[param] = np.random.normal(_mu, _sigma, num)
                if scale=='log':
                    _internalMC[param] = 10**_internalMC[param]

            _internalMC['PyriteSA'] = np.random.uniform(low=5.65, high=28.3, size=num)

            GP = GasProduction(_internalMC['bulk density'], _internalMC['Water-rock mass ratio'],
              _internalMC['U'], _internalMC['Th'], _internalMC['K'], _internalMC['PyriteSA'], _internalMC['Pyrite'],
              density_unit='ccm',
              K_unit = '%',
              sulfide_massfrac_unit='g/g',
              Energies=None, StoppingPower=None, Gs=None)

            GP_out = GP.output_df()
            # this is a df of the intermal MC of iteration i
            # we want to fit a normal to this, in log space

            for gas in out_gases:

                param_values = np.array(GP_out[gas])

                # remove any entries which went to -np.inf
                param_values = param_values[param_values > 0.]
                param_values = np.log10(param_values)

                density = spys.norm.fit(param_values)
                out_gases_dict[gas+' nom'].append(density[0])
                out_gases_dict[gas+' sigma'].append(density[1])

        # now out_gases_dict contains num entries, each of which is the
        # means and std devs of a fit in log space of the output gases.

        out_gases_df = pd.DataFrame(out_gases_dict)
        out_gases_df.to_csv('HardRockMC/data/DataOutput/CIS/'+fn)

        return out_gases_df
