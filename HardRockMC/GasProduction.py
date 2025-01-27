from copy import deepcopy
import pandas as pd
import ExceltoPandasMethods as EPM

class GasProduction:

    N_A = 6.02214076e23
    q = 1.60218e-19
    default_Energies = {
      'U' : {'a': 0.218e-3 , 'b': 0.146e-3, 'g': 0.113e-3},
      'Th' : {'a' : 0.061e-3, 'b': 0.027e-3, 'g': 0.048e-3},
      'K' : {'a': 0, 'b': 0.782e-7 , 'g': 0.243e-7}
    } # all in Gy /yr / ppm a.k.a. J/(kg yr ppm)

    default_StoppingPower = {'a': 1.5, 'b' : 1.25, 'g' : 1.14}


    default_G = {
      'H2':{'a':1.32, 'g' : 0.25, 'b':0.6}, # unit: molecules per 100 eV
      'SO4':{'a':2.1e-9, 'b':2.1e-9, 'g':2.1e-9} # unit: moles per m2 of sulfides per Gy
      }
    #H2 a, g from Sauvage et al 2021
    # b from Lin et al 2005 refs therein
    # SO4 same value for a, b, g, all from Li et al 2016, via Lefticariu et al 2010.

    @staticmethod
    def calc_W(bulk_density, porosity, water_density=1.):
        W_bulk = porosity*water_density / bulk_density
        W = W_bulk / (1-W_bulk)
        return W

    def __init__(self, bulk_density, W, U_ppm, Th_ppm, K, sulfide_SA, sulfide_massfrac,
      density_unit='', K_unit='', sulfide_massfrac_unit='',
      Energies=None, StoppingPower=None, Gs=None):


        MCdf =pd.DataFrame({ 'W': W,
          'U_ppm':U_ppm,
          'Th_ppm':Th_ppm,
          'sulfide_SA':sulfide_SA})
        if density_unit=='ccm':
            MCdf['bulk_density'] = bulk_density*1000
        elif density_unit=='kg/m3':
            MCdf['bulk_density'] = bulk_density
        else:
            raise ValueError('Please provide unit of bulk density. Options are ccm or kg/m3')

        if K_unit=='%':
            MCdf['K_ppm'] = K*1e4
        elif K_unit=='ppt':
            MCdf['K_ppm'] = K*1e3
        elif K_unit=='ppm':
            MCdf['K_ppm'] = K
        else:
            raise ValueError('Please provide unit of K concentration. Options are %, ppt, or ppm')

        if sulfide_massfrac_unit=='%':
            MCdf['sulfide_massfrac'] = sulfide_massfrac/100
        elif sulfide_massfrac_unit=='ppt':
            MCdf['sulfide_massfrac'] = sulfide_massfrac*1e-3
        elif sulfide_massfrac_unit=='ppm':
            MCdf['sulfide_massfrac'] = sulfide_massfrac*1e-6
        elif sulfide_massfrac_unit=='g/g':
            MCdf['sulfide_massfrac'] = sulfide_massfrac
        else:
            raise ValueError('Please provide unit of sulfide mass fraction. Options are %, ppt, or ppm')


        self.MCdf = MCdf

        if Energies == None:
            self.Energies = GasProduction.default_Energies
        else:
            self.Energies = Energies

        if StoppingPower == None:
            self.StoppingPower = GasProduction.default_StoppingPower
        else:
            self.StoppingPower = StoppingPower

        if Gs == None:
            self.G = GasProduction.default_G
        else:
            self.G = Gs


    def add_EnergyProductions_J_kg_yr(self):

        # total energy from alpha production; J/ (kg yr)
        self.MCdf['C_X__E_Xa'] = ((self.MCdf['U_ppm'] * self.Energies['U']['a']) +
          (self.MCdf['Th_ppm'] * self.Energies['Th']['a']) +
          (self.MCdf['K_ppm'] * self.Energies['K']['a']))

        # total energy from beta production; J/ (kg yr)
        self.MCdf['C_X__E_Xb'] = ((self.MCdf['U_ppm'] * self.Energies['U']['b']) +
          (self.MCdf['Th_ppm'] * self.Energies['Th']['b']) +
          (self.MCdf['K_ppm'] * self.Energies['K']['b']))

        # total energy from gamma production; J/ (kg yr)
        self.MCdf['C_X__E_Xg'] = ((self.MCdf['U_ppm'] * self.Energies['U']['g']) +
          (self.MCdf['Th_ppm'] * self.Energies['Th']['g']) +
          (self.MCdf['K_ppm'] * self.Energies['K']['g']))


    def H2Yield(self):

        """
        Compute He yield in moles / (m3 yr) using expressions in Warr et al 2023 Frontiers (and refs therein) based on Lin et al. (2005a&b).
        """

        H2_MCdf = deepcopy(self.MCdf) # for H2-specific working

        # molecules J / (bulk kg yr eV) from alpha
        H2_MCdf['isum_a'] = (H2_MCdf['C_X__E_Xa']*(H2_MCdf['W'] *self.StoppingPower['a']*self.G['H2']['a']) /
          (100*(1+(H2_MCdf['W'] *self.StoppingPower['a']))))

        # molecules J / (bulk kg yr eV) from beta
        H2_MCdf['isum_b'] = (H2_MCdf['C_X__E_Xb']*(H2_MCdf['W'] *self.StoppingPower['b']*self.G['H2']['b']) /
          (100*(1+(H2_MCdf['W'] *self.StoppingPower['b']))))

        # molecules J / (bulk kg yr eV) from gamma
        H2_MCdf['isum_g'] = (H2_MCdf['C_X__E_Xg']*(H2_MCdf['W'] *self.StoppingPower['g']*self.G['H2']['g']) /
          (100*(1+(H2_MCdf['W'] *self.StoppingPower['g']))))

        # convert to total H2 prod in moles per m3 per yr
        self.MCdf['Y_H2'] = (H2_MCdf['bulk_density'] /(GasProduction.q*GasProduction.N_A)) * (H2_MCdf['isum_a'] + H2_MCdf['isum_b'] + H2_MCdf['isum_g'])



    def SO4Yield(self):

        """
        Compute He yield in moles / (m3 yr) using expressions in Warr et al 2023 Frontiers (and refs therein) based on Lin et al. (2005a&b).
        """

        SO4_MCdf = deepcopy(self.MCdf) # for SO4-specific working

        # molecules J / (bulk kg yr eV) from alpha
        SO4_MCdf['isum_a'] = (SO4_MCdf['C_X__E_Xa']*(SO4_MCdf['W'] *self.StoppingPower['a']*self.G['SO4']['a']) /
          ((1+(SO4_MCdf['W'] *self.StoppingPower['a']))))

        # molecules J / (bulk kg yr eV) from beta
        SO4_MCdf['isum_b'] = (SO4_MCdf['C_X__E_Xb']*(SO4_MCdf['W'] *self.StoppingPower['b']*self.G['SO4']['b']) /
          ((1+(SO4_MCdf['W'] *self.StoppingPower['b']))))

        # molecules J / (bulk kg yr eV) from gamma
        SO4_MCdf['isum_g'] = (SO4_MCdf['C_X__E_Xg']*(SO4_MCdf['W'] *self.StoppingPower['g']*self.G['SO4']['g']) /
          ((1+(SO4_MCdf['W'] *self.StoppingPower['g']))))

        # convert to total SO4 prod in moles per bulk m3 per yr
        self.MCdf['Y_SO4'] =  (SO4_MCdf['isum_a'] + SO4_MCdf['isum_b'] + SO4_MCdf['isum_g']) * (SO4_MCdf['bulk_density']*SO4_MCdf['sulfide_SA']*SO4_MCdf['sulfide_massfrac'])


    def HeYield(self):
        """
        Compute He yield in moles / (m3 yr) using expressions in Warr et al 2023 Frontiers, based on Ballentine and Burnard 2002.
        """

        He_MCdf = deepcopy(self.MCdf)

        He_MCdf['U_contrib'] = (3.115e6+1.272e5)*He_MCdf['U_ppm']
        He_MCdf['Th_contrib'] = (7.710e5)*He_MCdf['Th_ppm']

        self.MCdf['Y_4He'] =  (He_MCdf['U_contrib'] + He_MCdf['Th_contrib']) * He_MCdf['bulk_density'] * 1e3 / GasProduction.N_A



    def ArYield(self):
        """
        Compute Ar yield in moles / (m3 yr) using expressions in Warr et al 2023 Frontiers, based on Ballentine and Burnard 2002.
        """

        self.MCdf['Y_40Ar'] = 102.2*self.MCdf['K_ppm'] * self.MCdf['bulk_density'] * 1e3 / GasProduction.N_A


    def output_df(self):


        self.add_EnergyProductions_J_kg_yr()
        self.H2Yield()
        self.SO4Yield()
        self.HeYield()
        self.ArYield()

        return self.MCdf


    # def get_uncertainties_df(self, save = False, preamble='data/DataOutput'):
    #
    #     """ possibly obsolete """
    #     split_columns = []
    #     for p in self.MCdf.columns.values.tolist():
    #         split_columns.append(p+' nom')
    #         split_columns.append(p+' sigma')
    #
    #     for p in self.MCdf.columns.values.tolist():
    #         EPM.create_nom_sigma(p, self.MCdf)
    #
    #     print(self.MCdf)
    #     unc_df = self.MCdf[split_columns]
    #     if type(save) == type(' '):
    #         unc_df.to_csv(preamble+save)
    #     return unc_df


    def save_output_df(self, fn, preamble='data/DataOutput/'):

        self.MCdf.to_csv(preamble+fn)
