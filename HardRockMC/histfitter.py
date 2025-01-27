
import matplotlib.pyplot as plt
import scipy

from scipy import stats as spys
from itertools import chain
import numpy as np
import math

from RetrieveProductionRates import RetrieveProductionRates
import ExceltoPandasMethods as EPM

class histfitter():

    """
    Class for generating propability distributions of parameters
    available in the ignace data sets
    """

    df_ref = {'K':'KThU_df',
      'bulk density' : 'D_df',
      'Th' : 'KThU_df',
      'Water-rock mass ratio' : 'W_df',
      'U' : 'KThU_df',
      'Pyrite' : 'KThU_df',
      'neighbor porosity' : 'W_df',
      'neighbor connected bulk mass [g/piece]' :'W_df',
      'surface area [scm/piece]' :'W_df',
      'volume [ccm/piece]' :'W_df',
      'Fracture porosity' : 'Bulk_df',
      'Y_H2' : 'output_df',
      'Y_4He' : 'output_df',
      'Y_40Ar' : 'output_df',
      'Y_SO4' : 'output_df',
      'Y_H2:Y_He' : 'output_df'}

    def __init__(self, BH, only_expt=None, exclude_expts=[],
      outputs = None,
      log_unit_adjust=0):
        """
        BH : str BH01, BH02, BH03 or all

        """

        self.KThU_df, self.D_df, self.W_df = EPM.get_dataframes(
          only_expt=only_expt, exclude_expts=exclude_expts)

        self.Bulk_df = EPM.get_bulk_dataframes()

        if BH != 'all':
            self.KThU_df = self.KThU_df.loc[self.KThU_df['Borehole']==BH]
            self.D_df = self.D_df.loc[self.D_df['Borehole']==BH]
            self.W_df  = self.W_df.loc[self.W_df['Borehole']==BH]
            self.Bulk_df  = self.Bulk_df.loc[self.Bulk_df['Borehole']==BH]

        self.output_df = None
        if outputs is not None:
            self.output_df = RetrieveProductionRates.ProductionRate_df(outputs)

        self.BH = BH
        self.log_unit_adjust=log_unit_adjust

        # if AQonly:
        #     self.W_df  = self.W_df.loc[self.W_df['Type']=='AQ'] # filter for only AQ water content






    def which_df(self, param):
        return eval('self.'+ self.df_ref[param])

    @staticmethod
    def find_best_fitter(param):
        return


    def get_param_values(self, param, scale='lin', min_infs_converge='min', return_num_infs=False):

        """
        return a numpy array of parameter measurements.
        If scale is 'log', return the parameter to the log base 10.

        min_infs_converge determines how to deal with -inf values encountered
        when using a log scale:
        if 'min', replace them with the smallest non-inf value
        if 'remove' remove them completely
        if a float, replace them with that value.

        Pass return_num_infs as True to return both the array of measurements
        and the number of infs encountered.
        """

        _df = self.which_df(param)
        param_values = None

        try:
            # try it with nominal values
            param_values = np.array(_df[param+' nom'].tolist())
        except:
            # try instead as-is.
            param_values = np.array(_df[param].tolist())

        if scale=='log':
            param_values = np.log10(param_values)
        elif scale == 'lin':
            pass
        else:
            raise ValueError('Undefined scale!')

        num_infs = None

        if min_infs_converge=='min':
            min_non_inf = np.nanmin(param_values[param_values != -np.inf])
            param_values[param_values < min_non_inf] = min_non_inf
            num_infs = (param_values == min_non_inf).sum()
        elif min_infs_converge=='remove':
            min_non_inf = np.nanmin(param_values[param_values != -np.inf])
            num_infs = (param_values <= min_non_inf).sum()
            param_values = param_values[param_values >= min_non_inf]
        elif type(min_infs_converge) == type(0.0):
            param_values[param_values < min_infs_converge] = min_infs_converge
            num_infs = (param_values == min_infs_converge).sum()

        param_values += self.log_unit_adjust

        if return_num_infs:
            return param_values, num_infs
        else:
            return param_values

    def fit(self, param, fitter, fitter_kwargs={}, xlims=None, scale='lin', min_infs_converge='min', xdim=1000):
        """
        Perform a fit with the distribution fitter.
        Return the density function, as well as x and y vals for a best fit
        line between limits xlims, if no xlims are passed, return a line
        between the max and min of the parameter range
        """

        _df = self.which_df(param)

        param_values = self.get_param_values(param, scale=scale, min_infs_converge=min_infs_converge)

        xmin, xmax = param_values.min(), param_values.max()
        xs = np.linspace(xmin,xmax, num=xdim)

        if xlims != None:
            xs= np.linspace(xlims[0], xlims[1], num=xdim)

        density = None
        if fitter == 'KDE':
            density = spys.gaussian_kde(param_values, **fitter_kwargs)
            return density, param_values, xs, density(xs)

        else:
            dist = eval("spys." + fitter)
            density = dist.fit(param_values, **fitter_kwargs)

            return density, param_values, xs, dist.pdf(xs, *density)


    def characterize(self, param, scale='lin', min_infs_converge='min'):

        # read in the data as an np array
        param_values = self.get_param_values(param, scale=scale, min_infs_converge=min_infs_converge)
        boxprops = {
          'median' : np.median(param_values),
          'mean' : np.mean(param_values),
          'min' : np.min(param_values),
          'pc5' : np.percentile(param_values, 5),
          'pc25' : np.percentile(param_values, 25),
          'pc75' : np.percentile(param_values, 75),
          'pc95' : np.percentile(param_values, 95),
          'max' : np.max(param_values),
          'param' : param,
          'BH' : self.BH,
          'N' : len(param_values)
        }
        return boxprops


    def resample(self, num, param, fitter, fitter_kwargs={}, scale='lin'):

        density, param_values, xs, pdf = self.fit(
           param, fitter, fitter_kwargs=fitter_kwargs, scale=scale
        )

        dist=None
        try:
            dist = eval("spys." + fitter).rvs(*density, size=num)
        except:
            dist = density.resample(num)[0]

        return dist
