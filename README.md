# denovoKemp

This repository contains a computational method for designing novel enzymes using a precomputed theozyme. The workflow is flexible and can be applied to any reaction where a theozyme is available, allowing the design of stable, functionally competent proteins with a focus on catalytic efficiency.

## Overview

The denovoKemp workflow combines several computational techniques to generate enzyme designs with structural diversity, stability, and catalytic activity. The process involves designing initial protein backbones, optimizing their active sites, and selecting the most promising candidates for further refinement. The method is based on the following key steps:

1. **Combinatorial Backbone Generation:**  
   We generate thousands of protein backbones by combining fragments from homologous proteins. A repositiory for this can be found [here](https://github.com/Fleishman-Lab/AbDesign_for_enzymes). These backbones then undergo [PROSS](https://pross.weizmann.ac.il/step/pross-terms/) design algorithm to stabilize the designed native state.

2. **Active Site Optimization:**  
   After backbone generation, we position the precomputed KE theozyme into each designed structure using geometric matching. We then refine the entire active site using Rosetta's atomistic calculations.

3. **Design Filtering:**  
   We generate millions of designs, which are filtered using an objective function that balances key factors for functional enzyme design.

4. **Further Stabilization:**  
   The top-scoring designs undergo additional stabilization steps using [FuncLib](https://ablift.weizmann.ac.il/step/fl_terms/) and [pSUFER](https://psufer.weizmann.ac.il/step/energy-threshold/), ensuring the active site and protein core are optimized for stability. 

## Methodology

The repository provides xmls, flag files and examples of command lines for each step. Follow the provided scripts in the different subdirectories. Inital pdb files are in the working_pdb directory.

- **Step 1:** Combinatorial assembly and design of protein backbones using homologous fragments, followed by PROSS design calculations for stability.
- **Step 2:** Geometric matching of the KE theozyme to each backbone structure and optimization of the active site using Rosetta.
- **Step 3:** Filtering of millions of designs based on energy and desolvation criteria to prioritize the most promising candidates.
- **Step 4:** Stabilization of the active site and protein core using FuncLib and pSUFER to refine the final designs.

## Requirements

- Python 3.x
- Jupyter notebook
- Rosetta

## Citation
If you use this method in your research, please cite the following:

High-efficiency Kemp eliminases by complete computational design
Dina Listov, Eva Vos, Gyula Hoffka, Shlomo Yakir Hoch, Andrej Berg, Shelly Hamer-Rogotner, Orly Dym, Shina Caroline Lynn Kamerlin, Sarel J. Fleishman
bioRxiv 2025.01.04.631280; doi: https://doi.org/10.1101/2025.01.04.631280
