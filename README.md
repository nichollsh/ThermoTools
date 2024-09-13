## Thermodynamics tools 

Includes tools for:
* Downloading and parsing data from the NIST-JANAF thermochemical 
tables.
* Calculating water phase change properties (Psat and enthalpy) using IAPWS95 
* Calculating saturation pressures for many gases using the Antoine equation
* Calculating enthalpy of phase change for many gases
* Compilation of thermodynamic properties into a single file (per gas)

### Setup

Simply run `pip install -e .` in your terminal, then use the notebooks.

### Sources

* NIST Webbook - https://webbook.nist.gov/
* JANAF Tables (Chase et al., 1985) - https://janaf.nist.gov/
* IAPWS95 / Wagner & Pruß (2002) - https://pubs.aip.org/aip/jpr/article/31/2/387/241937/The-IAPWS-Formulation-1995-for-the-Thermodynamic
* Feistel & Wagner (2006) - https://www.sciencedirect.com/science/article/pii/S0016703706020461
* Coker (2007) appendix C - https://doi.org/10.1016/B978-0-7506-7766-0.X5000-3
* Others documented in code

Copyright 2024 Harrison Nicholls under The Clear BSD License.
