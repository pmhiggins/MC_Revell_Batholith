from HardRockMC.MC_Implementor import MC_Implementor

number = 100000 # change to use a smaller sample size for the MC
date = '2024-09-18' # used as a prefix for output filenames


MCI = MC_Implementor()

"""
The MCI.implement function used below will run the MC model and save results
in the data/DataOutput directory.
"""

# df = MCI.implement(number, 'KDE', date+'_KDE.csv', 'nonweighted', porosity='primary')
# print(df)

# df = MCI.implement(number, 'KDE', date+'_KDE.csv', 'weighted', porosity='primary')
# print(df)

# df = MCI.implement(number, 'norm', date+'_norm.csv', 'nonweighted', porosity='primary')
# print(df)

# df = MCI.implement(number, 'norm', date+'_norm.csv', 'weighted', porosity='primary')
# print(df)


# Note: CIS will take a long time. For the results in the manuscript, 100000 iterations were used
# To make it shorter, change the `internal_num` kwarg for MC_Implementor.implement_CIS()
# df = MCI.implement(number, 'norm', date+'_norm.csv', 'CIS', porosity='primary')
# print(df)

"""
It is also possible to use secondary porosity, but it should be scaled
in log-space.
"""
# MCI.scales[0]='log'
# MCI.params[0]='Fracture porosity'
#
# df = MCI.implement(number, 'KDE', date+'_KDE_2ry.csv', 'nonweighted', porosity='secondary')
# print(df)
# df = MCI.implement(number, 'KDE', date+'_KDE_2ry.csv', 'weighted', porosity='secondary')
# print(df)
