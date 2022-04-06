# cliff_delineatools

## Tools for use with [CliffDelineaTool](https://github.com/zswirad/CliffDelineaTool) by [zswirad](https://github.com/zswirad).

Two short scripts for use in conjunction with [zswirad's CliffDelineaTool](https://github.com/zswirad/CliffDelineaTool), as described in Swirad Z.M. & Young A.P. 2021. Automating coastal cliff erosion measurements from large-area LiDAR datasets in California, USA. Geomorphology 389: 107799 (https://doi.org/10.1016/j.geomorph.2021.107799) for which the algorithm was originally developed (GitHub Zenodo link below).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5510191.svg)](https://doi.org/10.5281/zenodo.5510191)

### 1. points_from_dem.py
- Call the write_points() function to create points from a Digital Elevation Model (DEM) for use 
with cliff_delineatool.py 
- First creates transect lines across the DEM and then points along transect lines, pulling elevation information from
the DEM. 
- Writes Points to same directory as the DEM is in.

### 2. cliff_delineatool.py
- Adaptation of Zuzanna Swirad's CliffDelineaTool.
- Use the delineate_cliff() function to find the cliff base (and optionally top as well) for a
csv of cliff points created by points_from_dem.py module. 
- Writes both to the same folder as 
the cliff points csv.
- Variety of optional input args to change how cliff base is delineated.

### 3. polys_from_cliffpts.py
Call the build_polys() function to use cliff base points from cliff_delineatool.py to construct 
beach- and cliff-size MOP area tiles. 
 - First generates a line from the cliff base points.
 - Then uses the extent of the line to decide how many polygons to create.
 - Generates large polygons (of both beach + cliff) then splits beach/cliff using 
    cliff base line.
 - Assigns class code to beach/cliff polys: Beach = 1, Cliff = 0
 - Writes polygons and their attribute table to folder location of cliff base points input. 
