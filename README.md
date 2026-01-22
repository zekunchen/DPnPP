Code & Reproducibility Guide
This repository is provided for research purposes only. It includes the implementation of the UNNP framework and the core logic of the proposed DPnPP algorithm.

1. Quick Start (Reconstruction Demo)Direct Execution: Run Demo_Case1_UNNP.ipynb for an immediate example of image reconstruction using PyTorch.
* No License Required: The Jacobian matrix $J$ and necessary forward model data are pre-calculated via COMSOL and provided in the data/ folder3. Users can execute the reconstruction demo immediately without a COMSOL license.

2. Forward Problem & Jacobian Calculation
* The Jac/ folder contains the COMSOL + MATLAB simulation toolchain.

These files generate electric field data and reference conductivity. Users can recalculate the Jacobian matrix based on the formulas provided in the paper (please update the file paths in the scripts before use).

3. Algorithm Details
UNNP & DPnPP: This package provides the base Python implementation of UNNP. The core gradient calculation logic for DPnPP is illustrated in DPnPP_code.png. * Extended Reproducibility: Other benchmark algorithms mentioned in our paper (deepEIT, PINN, MLP, etc.) can be replicated by integrating the provided base code with the formulas in the manuscript  or by referring to the PnP-DIP source code(A Plug-and-Play Deep Image Prior https://github.com/zhaodongsun/pnp_dip/tree/master).

Citation
If this code or the DPnPP framework assists your research, please cite our paper:

"Deep Plug-and-play Prior for Enhanced Electrical Impedance Tomography"
