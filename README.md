
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

Describe the contents of the repository, and which pieces do what. You can use code text to refer to specific files or directories, like this: `a_file.ipynb`, `a_folder/`. In this template, the `figures.ipynb` contains a template Jupyter notebook for reproducing the figures accompanying the paper. Generally, the generated figures should be placed in the `figures/` directory.

If the analysis uses data that is maintained by a third party or stored separately (e.g., on Zenodo), then it can be linked to here. If the paper develops a method, then a test script should be included that implements the method, allowing users to verify the results before applying the method to their own data.

In general, local data should be organized in the `data/` directory and appropriate sub-folders. This should include raw data that we use (including from simulations) and processed data, which is saved separately.

Significant code should be in the `src/` directory. Code that is just used for making figures or minor data processing could be included here in the top-level directory.

Paper drafts, including cover letters, etc., can be placed in the `drafts/` directory. This directory also contains a `.gitignore` file that will make it such that the contents of this directory are not synced with GitHub.

To sync the above folders on GitHub, placeholder files have been placed in the `data/`, `figures/`, and `src/` directories. These files should be deleted when the template is copied for a project and files appear in each folder.

### Software dependencies

The code for generating and processing data is written in Julia (https://julialang.org/) and requires version 1.10 or later. Figures were created using Python 3.

# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
