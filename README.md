# adsorption
A python package for estimating the best isotherm (i.e., Langmuir, Freundlich, Temkin, Toth, Sips, DR) and kinetic models (i.e., PFO, PSO, Weber-Morris, Avrami, Bangham, Elovich) as well as Arrhenius parameters (i.e., Ea, A).

[Install](#Install) / [Usage](#Usage)

# Install
First, make sure you have a Python 3 environment installed.
To install from github:
```Python
pip install -e git+https://github.com/amvro23/adsorption/#egg=adsorption
```
Note: It might be useful to write "git+https://github.com/amvro23/adsorption/#egg=adsorption" if installing directly from a Python interpreter as # can be interpreted as a comment.

# Usage

```Python
from adsorption import (Isotherms, Kinetics, ModifiedArrhenius)
import numpy as np
import pandas as pd
```
Create an object with two required parameters x and y, where x represents the dimensionless equilibrium concentration [%] and y represents the equilibrium capacity [mg/g]. Default optional values are P = 1 atm, Mr = 44.01 g/mol for CO2 adsorption, T = 298.15 K, R = 8.205e-5 atm.m3/mol/K).
```Python
x = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])   
y = np.array([52.64, 72.59, 99.67, 141.79, 182.48, 203.68, 203.56, 204.33, 204.90])

isotherm = Isotherms(x, y)
```
You can obtain either a single model plot (e.g., Langmuir model).
```Python
isotherm.plot_langmuir_fit()
```
![langmuir](https://user-images.githubusercontent.com/91277572/208971022-c8579a81-6659-48c3-886b-a8dd92372c69.png)



