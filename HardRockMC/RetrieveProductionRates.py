import os, sys
import pandas as pd
import numpy as np


class RetrieveProductionRates:

    def ProductionRate_df(fn):

        loadas = '/data/DataOutput/'+fn

        fluxes_df = pd.read_csv(os.path.dirname(__file__)+loadas)

        try:

            fluxes_df['Y_H2:Y_He'] = fluxes_df['Y_H2']/fluxes_df['Y_4He']

        except:
            fluxes_df['Y_H2:Y_He nom'] = fluxes_df['Y_H2 nom']/fluxes_df['Y_4He nom']
            fluxes_df['Y_H2:Y_He sigma'] = fluxes_df['Y_H2 nom']*0.

        return fluxes_df


    def get_percentile(fn, pc, gas='H2'):

        _df = RetrieveProductionRates.ProductionRate_df(fn)

        p = np.array(fluxes_df['Y_'+gas+' nom'])
        print(fn, ' '+gas+' top +'(100-pc)+'% ', np.percentile(p, pc))

        return np.percentile(p, pc)
