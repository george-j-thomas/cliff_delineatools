# cliff_delineatools

## Tools for use with [CliffDelineaTool](https://github.com/zswirad/CliffDelineaTool) by [zswirad](https://github.com/zswirad).

Two short scripts for use in conjunction with [zswirad's CliffDelineaTool](https://github.com/zswirad/CliffDelineaTool), as described in Swirad Z.M. & Young A.P. 2021. Automating coastal cliff erosion measurements from large-area LiDAR datasets in California, USA. Geomorphology 389: 107799 (https://doi.org/10.1016/j.geomorph.2021.107799) for which the algorithm was originally developed (GitHub Zenodo link below).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5510191.svg)](https://doi.org/10.5281/zenodo.5510191)

### points_from_dem.py

- Call the write_points() function to create points from a Digital Elevation Model (DEM) for use 
with cliff_delineatool.py 
- First creates transect lines across the DEM and then points along transect lines, pulling elevation information from
the DEM. 
- Writes Points to same directory as the DEM is in.
