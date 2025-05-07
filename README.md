
This directory contains the code used to produce and analyse the results in the manuscript titled: *"Natural H2 and sulfate production via radiolysis in low porosity and permeability crystalline rocks"* by Higgins, Song, Warr and Sherwood Lollar (2025).

This code is able to replicate the Monte Carlo model presented, statistical analyses performed, and figures generated in that work and its supplemental material.

The `HardRockMC` directory contains code required to read in input parameters compiled from the NWMO reports (see references in the article, and the spreadsheet in `HardRockMC/data`), and perform a monte carlo simulation to estimate gas production of H2, SO4, Ar and He per unit of bulk rock (see article for methodology).

To facilitate replicability, the `MCmain.py` file exactly replicates the simulations performed in Higgins, Song et al. (2025), and the `Figures.py` file will regenerate all figures and supplement and save them to the Figures directory.
