
# Overview

This repository contains and data and scripts for reproducing the results accompanying the manuscript  

### Chronic infections can generate SARS-CoV-2-like bursts of viral evolution without epistasis
Edwin Rodr\'iguez-Horta<sup>1,2</sup>, John Strahan<sup>3</sup>, Aaron R. Dinner<sup>3</sup> and John P. Barton<sup>1,#</sup>

<sup>1</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine  
<sup>2</sup> Group of Complex Systems and Statistical Physics, Department of Theoretical Physics, Physics Faculty, University of Havana.  
<sup>3</sup> Department of Chemistry and James Franck Institute, University of Chicago  

<sup>#</sup> correspondence to [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)  

The preprint is available at __INSERT LINK HERE__.

# Contents

This repository contains code for evolutionary model simulations and a method to detect mutation bursts. The code is implemented in Julia and can be found in the `src/` directory. A demonstration of the code usage is available in `notebooks/Run_Test.ipynb`. 

Simulation results for each figure in the paper are located in the `data/` directory. Additionally, the `notebooks/` folder contains a collection of notebooks used to generate the figures for the paper.

### Software dependencies

The code for generating and processing data is written in Julia (https://julialang.org/) and requires version 1.10 or later. Figures were created using Python 3.

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
