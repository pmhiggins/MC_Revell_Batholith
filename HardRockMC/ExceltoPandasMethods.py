import os, sys; sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import pandas as pd
import math
import numpy as np
from uncertainties import ufloat as uf
from uncertainties import unumpy as unp
from uncertainties.unumpy import uarray as ua

"""
Collection of methods for translating the csv files in data/DataInput into a pandas dataframe,
which can be used for the calculations in python. Can also send output to a csv file (data/DataProcessed)
for easier retrieval / analysis in other softwares such as Excel.

"""

def from_csv(filename):
    df = pd.read_csv(filename)
    return df

def create_nom_sigma(param, _df):
    """
    create new columns in database _df corresponding to the nominal and
    1 sigma uncertainties in parameter param
    """
    _df[param+' nom'], _df[param+' sigma'] = unp.nominal_values(_df[param]), unp.std_devs(_df[param])


def calc_dens_unc(D_df):
    """
    Read in data from the density dataframe D_df, create extra columns for
    uncertainties, fill in gaps where possible, create columns of
    uncertainties.ufloat objects for error propagation
    Return newly populated dataframe

    Data in D_df should contain:
        AQ, PW (from early reports),
        BD,TP (BH01 interim report),
        PP (Alberta reports)

    """

    # add ufloat columns of those where we have nominals and sigmas in input
    D_df['bulk mass [g]'] = unp.uarray(D_df['bulk mass b.e. [g] nom'], D_df['bulk mass b.e. [g] sigma'])
    D_df['dry mass [g]'] = unp.uarray(D_df['dry mass [g] nom'], D_df['dry mass [g] sigma'])
    D_df['diameter core [cm]'] = unp.uarray(D_df['diameter core [cm] nom'], D_df['diameter core [cm] sigma'])
    D_df['height core [cm]'] = unp.uarray(D_df['height core [cm] nom'], D_df['height core [cm] sigma'])

    # get radius column
    D_df['radius core [cm]'] = D_df['diameter core [cm]'] / 2.

    # for BD and TP, create ufloat entry directly from volumes
    # we do this to calculate the apparent 'height' of the core piece,
    # because that information is not available
    D_df.loc[(D_df['Type']=='BD') | (D_df['Type']=='TP'), 'volume core [ccm]'] = unp.uarray(
      D_df.loc[(D_df['Type']=='BD') | (D_df['Type']=='TP'), 'volume core [ccm] nom'],
      D_df.loc[(D_df['Type']=='BD') | (D_df['Type']=='TP'), 'volume core [ccm] sigma'])

    # now set their height
    D_df.loc[(D_df['Type']=='BD') | (D_df['Type']=='TP'), 'height core [cm]'] = D_df['volume core [ccm]'] / (math.pi*D_df['radius core [cm]']*D_df['radius core [cm]'])
    create_nom_sigma('height core [cm]', D_df)

    # now volume ufloat can be created for all entries, because they all have height data
    D_df['volume core [ccm]'] = D_df['radius core [cm]']*D_df['radius core [cm]']*D_df['height core [cm]']*math.pi
    create_nom_sigma('volume core [ccm]', D_df)

    # likewise for surface area
    D_df['surface area core [scm]'] = (2*math.pi*D_df['radius core [cm]']*D_df['radius core [cm]']) + (D_df['height core [cm]']*math.pi*2*D_df['radius core [cm]'])
    create_nom_sigma('surface area core [scm]', D_df)

    # and finally bulk and dry density
    D_df['bulk density'] = D_df['bulk mass [g]'] / D_df['volume core [ccm]']
    create_nom_sigma('bulk density', D_df)

    D_df['dry density'] = D_df['dry mass [g]'] / (D_df['volume core [ccm]'])
    create_nom_sigma('dry density', D_df)

    # obsolete code commented below, we used to have a porosity column here
    # but that gets confusing with the W_df dataframe

    # D_df['calc. phi'] = (D_df['m dry surf b.e [g]'] - D_df['m dry [g]']) / (D_df['radius core [cm]']*D_df['radius core [cm]']*D_df['height core [cm]']*math.pi*uf(1,0.03))
    # D_df['calc. phi nom'], D_df['calc. phi sigma'] = unp.nominal_values(D_df['calc. phi']), unp.std_devs(D_df['calc. phi'])

    return D_df

def add_KThU_err(df, pc_err=10.):
    """ add an arbitraty % error pc_error to the concentration of radielements """
    df['K sigma'] = df['K nom']*pc_err/100
    df['U sigma'] = df['U nom']*pc_err/100
    df['Th sigma'] = df['Th nom']*pc_err/100
    return df


def collate_densities(df_from, df_to):

    """
    possibly obsolete, used to make sure the densities between one and another
    dataframe were the same
    """

    samples = df_from['Sample']
    df_to["calc. wet density nom"] = np.nan
    df_to["calc. wet density sigma"] = np.nan

    df_to["calc. dry density nom"] = np.nan
    df_to["calc. dry density sigma"] = np.nan

    for s in samples:
        for d in ['calc. wet density', 'calc. dry density', 'volume core [ccm]', 'surface area core [scm]']:
            _dens = df_from[(df_from['Sample'] == s)][d].tolist()[0]
            indicies = df_to.index[(df_to['Sample'] == s)]
            for i in indicies:
                df_to.loc[i,d+' nom'] = _dens.n
                df_to.loc[i,d+' sigma'] = _dens.s


    df_to["calc. wet density"] = unp.uarray(df_to["calc. wet density nom"], df_to["calc. wet density sigma"] )
    df_to["calc. dry density"] = unp.uarray(df_to["calc. dry density nom"], df_to["calc. dry density sigma"] )
    df_to["volume core [ccm]"] = unp.uarray(df_to["volume core [ccm] nom"], df_to["volume core [ccm] sigma"] )

    return df_from, df_to

def setup_porosity_df(W_df, D_df):

    """
    Using the density dataframe and water-loss measurements, estimate the
    effective porosity from the Ignace experiments, even for those where we
    do NOT have volume (therefore density) estimates by approximating the density
    as either that of the closest "neighbor" or using a mean value for the whole dataset.
    """

    W_df['bulk mass [g]'] = unp.uarray(W_df['bulk mass [g] nom'], W_df['bulk mass [g] sigma'])
    create_nom_sigma('bulk mass [g]', W_df)

    # where available in the input (everything except PW, create ufloats for volume
    W_df['volume [ccm/piece]'] = unp.uarray(W_df['volume [ccm/piece] nom'], W_df['volume [ccm/piece] sigma'])
    # where available in the input (PWoutdiff, BD, TP, PP), create ufloats for surface area
    W_df['surface area [scm/piece]'] = unp.uarray(W_df['surface area [scm/piece] nom'], W_df['surface area [scm/piece] sigma'])
    #  calc SA of those 'sphere' (e.g. not a cylinder core, but crushed pieces)
    # Currently: AQ, IsoEx
    W_df.loc[W_df['shape'] == 'sphere', 'surface area [scm/piece]'] = 4*math.pi*((3*W_df['volume [ccm/piece]']/(4*math.pi))**(2/3))

    # now we have volume and SA for all pieces, except PW.
    # currently PW is ignored - need to think about what to do here.
    # volume and SA not known, 'end pieces which then had 20mm shaved off' so COULD assume cylinder with avg radius?
    create_nom_sigma('surface area [scm/piece]', W_df)

    # uncertainty in W (mass connected water / (bulk or dry) mass rock) is tricky. first calculate dry mass:
    W_df['dry mass [g] nom'] = W_df['bulk mass [g] nom']*(1-W_df['W wet'])
    # erro est to 0.05 g here, from the reports
    W_df['dry mass [g] sigma'] = 0.05
    W_df['dry mass [g]'] = unp.uarray(W_df['dry mass [g] nom'], W_df['dry mass [g] sigma'])

    # now re-calculate W dry and W wet with the masses and their uncertainties
    W_df['W wet'] = (W_df['bulk mass [g]'] - W_df['dry mass [g]']) / W_df['bulk mass [g]']
    W_df['Water-rock mass ratio'] = (W_df['bulk mass [g]'] - W_df['dry mass [g]']) / W_df['dry mass [g]']

    create_nom_sigma('W wet', W_df)
    create_nom_sigma('Water-rock mass ratio', W_df)


    # for eventual use with the bonus BH 1 pieces?
    diam = uf(6.09492307692308, 0.0204727774223386)

    # time for 'neighbor' and 'general' densities
    samples = D_df['Sample']
    W_df["neighbor bulk density"] = np.full(len(W_df), np.nan, dtype=object)
    W_df["neighbor dry density"] = np.full(len(W_df), np.nan, dtype=object)
    W_df["BHmean dry density"] = np.full(len(W_df), np.nan, dtype=object)
    W_df["BHmean bulk density"] = np.full(len(W_df), np.nan, dtype=object)


    # look through the D_df dataframe and assign the neighbor denisties accordingly
    for s in samples:
        for d in ['bulk density', 'dry density']:
            # find the entry in D_df corresponding to the same sample
            # aka big core piece the indidual experimental pieces were taken from
            _dens = D_df[(D_df['Sample'] == s)][d].tolist()[0]
            # find the indicies in W_df it corresponds to
            indicies = W_df.index[(W_df['Sample'] == s)]
            # assign the density in W_df
            for i in indicies:
                W_df.loc[i,'neighbor '+d] = _dens

    # also make a 'mean' density per borehole, using all the calculations and
    # finding the std deviation to use as 1 sigma
    BHs = ['BH01', 'BH02', 'BH03']
    for BH in BHs:
        for d in ['bulk density', 'dry density']:
            mean_dens = np.array(D_df[(D_df['Borehole'] == BH)][d+' nom'].tolist()).mean()
            std_dens = np.array(D_df[(D_df['Borehole'] == BH)][d+' sigma'].tolist()).std()

            W_df.loc[W_df['Borehole'] == BH, 'BHmean '+d] = uf(mean_dens, std_dens)

    #now create the nom and std columns for these new densities
    create_nom_sigma("neighbor bulk density", W_df)
    create_nom_sigma("neighbor dry density", W_df)
    create_nom_sigma("BHmean bulk density", W_df)
    create_nom_sigma("BHmean dry density", W_df)


    # finally, each data point can now get its own 'neighbor porosity' and 'general porosity'
    # may be worth rethinking the uncertainty in water density? This was in one of the reports....
    # NB all effective porosities are estimated with the same expression:
    # phi = W_bulk * rho_bulk / rho_w
    W_df['neighbor porosity'] = W_df['neighbor bulk density'] * W_df['W wet'] /uf(1,0.03)
    create_nom_sigma('neighbor porosity', W_df)
    W_df['general porosity'] = W_df['BHmean bulk density'] * W_df['W wet'] / uf(1,0.03)
    create_nom_sigma('general porosity', W_df)


    # finally get the 'connected mass' i.e. the mass of each 'piece'
    # this is useful for comparisons in a similar way to the volume or surface area per piece
    # but it is based on density so also has neighbor and general possibilities!
    W_df['neighbor connected bulk mass [g/piece]'] = W_df['neighbor bulk density'] * W_df['volume [ccm/piece]']
    create_nom_sigma('neighbor connected bulk mass [g/piece]', W_df)

    W_df['general connected bulk mass [g/piece]'] = W_df['BHmean bulk density'] * W_df['volume [ccm/piece]']
    create_nom_sigma('general connected bulk mass [g/piece]', W_df)

    return W_df



def get_dataframes(only_expt=None, exclude_expts=[]):

    """
    Main interface of EPM with rest of code. Formulates the dataframes as
    needed and returns them.
    """
    KThU_df = from_csv(os.path.dirname(__file__)+'/data/DataInput/NWMO_data/KThU_BH.csv')
    DBH3_df = from_csv(os.path.dirname(__file__)+'/data/DataInput/NWMO_data/D_df.csv')
    PBH3_df = from_csv(os.path.dirname(__file__)+'/data/DataInput/NWMO_data/AllWaterContent.csv')

    D_df = calc_dens_unc(DBH3_df)
    KThU_df = add_KThU_err(KThU_df, pc_err=10.)
    W_df = setup_porosity_df(PBH3_df, D_df)

    if only_expt != None:
        W_df = W_df[W_df['Type'] == only_expt]
        D_df = D_df[D_df['Type'] == only_expt]

    for e in exclude_expts:
        W_df = W_df[W_df['Type'] != e]
        D_df = D_df[D_df['Type'] != e]


    # D_df, W_df = collate_densities(D_df, PBH3_df)
    # EPM.to_csv(D_df, 'data/BH03manual_DensityPorosity.csv')
    return KThU_df, D_df, W_df

def get_bulk_dataframes():
    Bulk_df = pd.read_csv(os.path.dirname(__file__)+'/data/DataInput/NWMO_data/FracturePorosity.csv')
    return Bulk_df

def take_logs(dfs, logparams):
    """ take logs of indicated columns in  passed list of dataframes """
    for p in logparams:
        for df in dfs:
            if p in df.columns:
                df[p] = np.log10(df[p])
                print('Taking logs of '+p)
            elif p+' nom' in df.columns:
                df[p+' nom'] = np.log10(df[p+' nom'])
                print('Taking logs of '+p)

    return dfs

def to_csv(df, filename):
    df.to_csv(filename)



def distributions_to_csv(names, distros, filename='data/distribution.csv', negpurge=True):
    _dict = {n:d for n, d in zip(names, distros)}
    df = pd.DataFrame(_dict)
    if negpurge:
        df[df < 0] = 0
    # df['porosity [Vw / Vb]'] = df['bulk density [g/ccm]']*df['W [m_w/m_b]']
    df.to_csv(filename)
    return df


# KThU_df, D_df, W_df = get_dataframes()
# D_df.to_csv('temp.csv')
