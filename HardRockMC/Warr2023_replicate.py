import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import ExceltoPandasMethods as EPM
from SampleGenerator import SampleGenerator
from GasProduction import GasProduction
import sigfig


"""
This is an example of using HardRockMC to estimate gas production somewhere
other than Revell, here replicating the results from Warr et al (2023)
for Kidd Creek.
"""

SG = SampleGenerator(100000, 'norm', params=[ 'neighbor porosity',
 'bulk density',
 'U',
 'Th',
 'K',
 'water density',
 'Pyrite',
 'PyriteSA'])


# generates equivalent distributions to warr et al 2023 for input into the MC and saves them
# to the filename fn in csv/DataInput/Redistributed
fn_norm = SG.GenerateWarr2023Samples(100000, type='norm')

SG = SampleGenerator(100000, 'norm', params=[ 'neighbor porosity',
 'bulk density',
 'U',
 'Th',
 'K',
 'water density',
 'Pyrite',
 'PyriteSA'])


fn_uniform = SG.GenerateWarr2023Samples(100000, type='uniform')


# now extract those distributions again
loadas_norm = os.path.dirname('data/DataInput/Redistributed')+'/'+fn_norm
loadas_uniform = os.path.dirname('data/DataInput/Redistributed')+'/'+fn_uniform

dist_df_norm = EPM.from_csv(loadas_norm)
dist_df_uniform = EPM.from_csv(loadas_uniform)


# Create a GasProduction object
GP_norm = GasProduction(dist_df_norm['bulk density'], dist_df_norm['Water-rock mass ratio'], dist_df_norm['U'], dist_df_norm['Th'], dist_df_norm['K'], dist_df_norm['PyriteSA'], dist_df_norm['Pyrite'],
  density_unit='ccm',
  K_unit = 'ppm',
  sulfide_massfrac_unit='%',
  Energies=None, StoppingPower=None, Gs=None)

GP_uniform = GasProduction(dist_df_uniform['bulk density'], dist_df_uniform['Water-rock mass ratio'], dist_df_uniform['U'], dist_df_uniform['Th'], dist_df_uniform['K'], dist_df_uniform['PyriteSA'], dist_df_uniform['Pyrite'],
  density_unit='ccm',
  K_unit = 'ppm',
  sulfide_massfrac_unit='%',
  Energies=None, StoppingPower=None, Gs=None)

# this calucates the gas production rates and adds them to its internal DataFrame
# Sulfate production is coded up, included, but not formally reviewed yete
GP_norm.output_df()
GP_uniform.output_df()

# save the output to a new csv of the same name,
# but now in directory csv/DataProcessed
# GP_norm.save_output_df(fn_norm)
# GP_uniform.save_output_df(fn_uniform)

for _GP, _type in zip([GP_norm, GP_uniform], ['normal', 'uniform']):
    print(' ####################   '+_type+'     ####################')
    for p in ['Y_SO4', 'Y_H2', 'Y_40Ar', 'Y_4He']:
        mean = sigfig.round(_GP.MCdf[p].mean(),3)
        std = sigfig.round(_GP.MCdf[p].std(),3)

        print(p, mean, std)
