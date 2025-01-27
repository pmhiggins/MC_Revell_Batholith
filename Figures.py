"""
This code will generate every figure in the manuscript Higgins, Song et al. (submitted).
"""
import matplotlib.pyplot as plt

""" Figure 2 (porosity) """
from Analysis.inputs.PorosityScale import PorosityScale
# PorosityScale.allBH('F2_Porosity_volume.pdf')


""" Figure 3, S2 and S3 """
from Analysis.inputs.BHbyBH_inputcomparison import BHbyBH_inputcomparison

# BHbyBH_inputcomparison.Boxplots(all_model=['KDE', 'norm'], BH_model='KDE',
#   params = ['neighbor porosity','Fracture porosity','bulk density','U','Th','K'],
#   filename='F3_input_distributions_boxplots.pdf', skew=True)

# BHbyBH_inputcomparison.Boxplots(all_model=['KDE', 'norm'], BH_model='norm',
#   params = ['neighbor porosity','Fracture porosity','bulk density','U','Th','K'],
#   filename='FS2_input_distributions_boxplots.pdf')
#
# BHbyBH_inputcomparison.Boxplots(all_model=['KDE', 'norm'], BH_model='KDE',
#   params = ['Pyrite', 'Water-rock mass ratio'],
#   filename='FS3a_input_distributions_boxplots.pdf')
#
# BHbyBH_inputcomparison.Boxplots(all_model=['KDE', 'norm'], BH_model='norm',
#   params = ['Pyrite', 'Water-rock mass ratio'],
#   filename='FS3b_input_distributions_boxplots.pdf')


""" KS tests (Table S2) """
# BHbyBH_inputcomparison.KS_test(
#   params = ['Water-rock mass ratio', 'neighbor porosity','Fracture porosity','bulk density','U','Th','K'])

""" general statistical properties for inputs (Table S1) """
# BHbyBH_inputcomparison.input_data_basic_properties(
#   params = ['Water-rock mass ratio', 'neighbor porosity','Fracture porosity','bulk density','U','Th','K', 'Pyrite'])


""" Figure S1 """
from Analysis.inputs.depth_scatter import depth_scatter
# depth_scatter.plot_depth_ind_joints(means=True, intervals=100)

""" Figure S4 """
from Analysis.inputs.density_vs_W import density_vs_W
# density_vs_W.dens_vs_W()


""" Figure S5 """
from Analysis.inputs.normal_uncertainties import normal_uncertainties
# normal_uncertainties.normal_confidencelimits(
#   params=['neighbor porosity', 'bulk density', 'U', 'Th', 'K', 'Pyrite'])


""" Figure 4 and Figure S6 """
date = '2024-09-18'
from Analysis.outputs.gas_flux import gas_flux

# # Figure 4
# gas_flux.individual_output_boxes(fn_KDE=date+'_KDE.csv',
#   fn_norm=date+'_norm.csv', fn_CIS=date+'_norm.csv',
#   params = {'Y_H2': r'H$_2$ production rate [mol (m$^{3}$ rock yr)$^{-1}$]',
#     'Y_SO4':  r'Sulfate production rate [mol (m$^{3}$ rock yr)$^{-1}$]'},
#   weighted=False,
#   CIS_xdim=1000,
#   fig_preamble='F4__',
#   scale='log')
#
# # Figure S6
# gas_flux.individual_output_boxes(fn_KDE=date+'_KDE.csv',
#   fn_norm=date+'_norm.csv', fn_CIS=date+'_norm.csv',
#   params = {'Y_H2': r'H$_2$ production rate [mol (m$^{3}$ rock yr)$^{-1}$]',
#     'Y_SO4':  r'Sulfate production rate [mol (m$^{3}$ rock yr)$^{-1}$]'},
#   weighted=True,
#   CIS_xdim=1000,
#   fig_preamble='FS6_',
#   scale='log')
#
""" Basic stats params (Supplemental Table S3) """
# gas_flux.output_data_basic_properties(
#   params = ['Y_H2', 'Y_SO4'],
#   fn_KDE=date+'_KDE.csv',
#   fn_norm=date+'_norm.csv')
