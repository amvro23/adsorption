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
You can obtain either a single model plot (e.g., Langmuir model, Freundlich model).
```Python
isotherm.plot_langmuir_fit()
```
![langmuir](https://user-images.githubusercontent.com/91277572/208971022-c8579a81-6659-48c3-886b-a8dd92372c69.png)
```Python
isotherm.plot_freundlich_fit()
```
![freundlich](https://user-images.githubusercontent.com/91277572/208971538-1b98051e-b8c1-47ad-b55d-8837146c31ed.png)
Or you can obtain all fitting models.
```Python
isotherm.plot_all_models()
```
![all_isotherms](https://user-images.githubusercontent.com/91277572/208971930-40142a78-459c-4e70-840b-88829d8ffe2a.png)

You can assess the fitting of isotherm models.
```Python
isotherm.assess_fit()
```
```Python
Out
{'Langmuir R2': 0.9853798413181968,
 'Freundlich R2': 0.7822175305642807,
 'Temkin R2': 0.8766497302509872,
 'Toth R2': 0.9852294550091175,
 'Sips R2': 0.9851843142530591,
 'DR R2': 0.8808791659298837}
```
 The package can directly return the best fitting model.
```Python
isotherm.best_fit()
```
```Python
Out
The best model is that of Langmuir R2 = 0.9853798413181968
```
You can also obtain a dataframe with all the calculated parameters with their units.
```Python
isotherm.all_params()
```
```Python
                                           Parameters               Values
0                                   K_Langmuir [L/mg]  0.05971152127312045
1                                qmax_Langmuir [mg/g]   211.83876099132755
2   K_Freundlich [L^(1/n_Freundlich)·mg^(1-1/n_Fre...     71.6650763022144
3                                        n_Freundlich    6.254943331258903
4                                     A_Temkin [L/mg]   2.9827960331000063
5                                    B_Temkin [J/mol]     90.1390756026659
6                             K_Toth [mg^z·L^-n_toth]   15.637639450829793
7                                              n_Toth   0.9832070920173336
8                                    qmax_Toth [mg/g]   212.27301879029434
9                        K_Sips [L^n_Sips·mg^-n_Sips]   0.0685096207219515
10                                             n_Sips    1.056671824783758
11                                   qmax_Sips [mg/g]    213.6471311521534
12                                       E_DR [J/mol]   147.67076438840127
13                                     qmax_DR [mg/g]    192.6633619518013
```

