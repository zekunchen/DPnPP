This code is provided for research purpose only.

The "Jac" folder contains COMSOL + MATLAB simulation files. These can generate electric field data and reference conductivity. Afterwards, you can calculate the Jacobian matrix based on the formula (please modify the path before use; the files are already complete). The Jacobian matrix J and forward model data are pre-calculated using COMSOL and provided in the data folder. Users can run the reconstruction directly without a COMSOL license.

Run "Demo_Case1_UNNP.ipynb" for an example. (Using the pytorch)

This code package is a Python implementation of UNNP. DPnPP is implemented based on this UNNP code (DPnPP_code.png is the core code for gradient calculation). Other algorithms such as deepEIT, PINN, MLP and others mentioned in the papers can all be replicated based on the provided basic code. Please reproduce it according to this code and the formulas in the paper. Or based on the source code provided in this paper, “A Plug-and-Play Deep Image Prior.(https://github.com/zhaodongsun/pnp_dip/tree/master)”

Please cite our paper if this code is used to motivate any publications.

"Deep Plug-and-play Prior for Enhanced Electrical Impedance Tomography" 
