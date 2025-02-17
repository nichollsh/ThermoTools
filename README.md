## Thermodynamics tools

Includes tools for:
* Downloading and parsing heat capacities from the JANAF tables
* Calculating water phase change properties (Psat and enthalpy) using IAPWS95
* Calculating saturation pressures for using the Antoine and Clausius-Clapeyron equations
* Tabulating enthalpy / latent heat of phase change
* Tabulating equations of state (Van der Waals, AQUA)
* Compilation of thermodynamic properties into a single file for each species

### Setup

Simply run `pip install -e .` in your terminal, then use the `workflow/` notebooks in order.

### Sources

* NIST Webbook - https://webbook.nist.gov/
* JANAF Tables (Chase et al., 1985) - https://janaf.nist.gov/
* IAPWS95, Wagner & Pruß (2002) - https://pubs.aip.org/aip/jpr/article/31/2/387/241937/The-IAPWS-Formulation-1995-for-the-Thermodynamic
* Feistel & Wagner (2006) - https://www.sciencedirect.com/science/article/pii/S0016703706020461
* Coker (2007), Appendix C - https://doi.org/10.1016/B978-0-7506-7766-0.X5000-3
* AQUA water equation of state - https://www.aanda.org/articles/aa/abs/2020/11/aa38367-20/aa38367-20.html 
* CMS19 hydrogen equation of state - https://iopscience.iop.org/article/10.3847/1538-4357/aaf99f
* Van der Waals coefficients - https://en.wikipedia.org/wiki/Van_der_Waals_constants_(data_page)
* Others documented in code

Copyright (c) 2025 Harrison Nicholls, under The Clear BSD License.
