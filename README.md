# adsorption
A python package for estimating the best isotherm (i.e., Langmuir, Freundlich, Temkin, Toth, Sips, DR) and kinetic models (i.e., PFO, PSO, Weber-Morris, Avrami, Bangham, Elovich) as well as Arrhenius parameters (i.e., Ea, A).

[Install](#Install) / [Usage](#Usage) /  [Isotherms](#Isotherms) / [Kinetics](#Kinetics) / [Arrhenius](#Arrhenius) / [References](#References)

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
## Isotherms

Create an object with two required parameters x and y, where x represents the dimensionless equilibrium concentration [%] and y represents the equilibrium capacity [mg/g]. Default optional values are P = 1 atm, Mr = 44.01 g/mol for CO2 adsorption, T = 298.15 K, R = 8.205e-5 atm.m3/mol/K).
```Python
x = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])   
y = np.array([52.64, 72.59, 99.67, 141.79, 182.48, 203.68, 203.56, 204.33, 204.90])

isotherm = Isotherms(x, y)
```
You can obtain either a single isotherm model plot (e.g., Langmuir model, Freundlich model).
```Python
isotherm.plot_langmuir_fit()
```
![langmuir](https://user-images.githubusercontent.com/91277572/208971022-c8579a81-6659-48c3-886b-a8dd92372c69.png)
```Python
isotherm.plot_freundlich_fit()
```
![freundlich](https://user-images.githubusercontent.com/91277572/208971538-1b98051e-b8c1-47ad-b55d-8837146c31ed.png)

Or you can obtain all isotherm models.
```Python
isotherm.plot_all_models()
```
![all_isotherms](https://user-images.githubusercontent.com/91277572/208971930-40142a78-459c-4e70-840b-88829d8ffe2a.png)

You can assess the isotherm models.
```Python
isotherm.assess_fit()
```
```
Out
{'Langmuir R2': 0.9853798413181968,
 'Freundlich R2': 0.7822175305642807,
 'Temkin R2': 0.8766497302509872,
 'Toth R2': 0.9852294550091175,
 'Sips R2': 0.9851843142530591,
 'DR R2': 0.8808791659298837}
```
 The package can directly return the best isotherm model.
```Python
isotherm.best_fit()
```
```
Out
The best model is that of Langmuir R2 = 0.9853798413181968
```
You can also obtain a dataframe with all the calculated parameters of isotherm equations with their units.
```Python
isotherm.all_params()
```
```
Out
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
Or you can have access in isotherm parameters individually in order to create your own plots.
```Python
isotherm.dubinin_radushkevich_params()
```
```
Out
{'E_DR [J/mol]': 147.67076438840127, 
'qmax_DR [mg/g]': 192.6633619518013}
```
```Python
isotherm.sips_curve(x)
```
```
Out
array([ 0.04084583,  0.09719205,  0.18721156,  0.36046663,  0.85595811,
        1.64334054,  3.14432104,  7.33485598, 13.69841098])
```
## Kinetics
Create an object with two required parameters x and y, where x represents the adsorption time [min] and y represents the cumulative adsorption capacity [mg/g]. In this example a csv file was used to obtain these values.
```Python
df = pd.read_csv('adsorption_kinetics.csv')
x = df.loc[:, 'minutes'].values
y = df.loc[:, 'qt'].values

kinetic = Kinetics(x, y)
```
You can obtain either a single kinetic model plot (e.g., Bangham model).
```Python
kinetic.plot_bangham_fit()
```
![bangham](https://user-images.githubusercontent.com/91277572/208979002-cfb89967-123a-498c-b9ff-423ceb0dd7c3.png)

Or you can obtain all kinetic models.
```Python
kinetic.plot_all_models()
```
![kinetic_models](https://user-images.githubusercontent.com/91277572/208979310-532873b0-2073-4a32-92f4-730991f1ff28.png)

You can assess the kinetic models.
```Python
kinetic.assess_fit()
```
```
Out
{'PFO R2': 0.9862768756329133,
 'PSO R2': 0.9787375821197556,
 'WEBER-MORRIS R2': 0.9568081287159736,
 'AVRAMI R2': 0.9967027963259552,
 'BANGHAM R2': 0.9983438488830458,
 'ELOVICH R2': nan}
 ```
 The package can directly return the best kinetic model.
 ```Python
kinetic.best_fit()
```
 ```
Out
 The best model is that of BANGHAM R2 = 0.9983438488830458
 ```
You can also obtain a dataframe with all the calculated parameters of kinetic equations with their units.
 ```Python
kinetic.all_params()
```
 ```
Out
                 Parameters                  Values
0             k_pfo [1/min]    0.018620160753212878
1           qmax_PFO [mg/g]       74.01352655075553
2    k_pso [g mg^-1 min^-1]  0.00012514880623542925
3           qmax_PSO [mg/g]      109.00232132690607
4   k_wm [mg g^-1 min^-0.5]       6.184044312110834
5                         C   6.836063566443158e-18
6          k_avrami [1/min]    0.035198182776525935
7        qmax_avrami [mg/g]       65.19974372878553
8                  n_avrami      1.7755718014035597
9       k_bangham [1/min^n]   0.0047733638932455965
10     qmax__bangham [mg/g]       63.29831524175902
11                n_bangham      1.4416628511124052
12     a_elovich [mg/g/min]       4.602398573502221
13        b__elovich [g/mg]     0.05477589170964753
 ```
 Or you can have access in kinetic parameters individually in order to create your own plots.
  ```Python
kinetic.avrami_params()
```
 ```
Out
{'k_avrami [1/min]': 0.035198182776525935,
 'qmax_avrami [mg/g]': 65.19974372878553,
 'n_avrami': 1.7755718014035597}
```
  ```Python
kinetic.pso_curve(x)
```
  ```
  Out
array([ 0.        ,  0.08073146,  0.16136812, ..., 67.24288195,
       67.25473634, 67.266584  ])
```
## Arrhenius
After finding the best kinetic model, the equation is applied to the data obtained at 3 to 4 different temperatures (at least) in order to create the following object, where x represents the temperatures in which adsorption tests were performed and y represents the k value of the best kinetic model at the same temperatures. In this example Bangham model was the best to describe the process of CO2 adsorption.
  ```Python
df_kin = pd.DataFrame({'Value':[0.00478, 0.00583, 0.00728, 0.01956],
                      'Temperature':[298.15, 308.15, 323.15, 373.15]})

x = df_kin.loc[:, 'Temperature'].values
y = df_kin.loc[:, 'Value'].values

arrh = ModifiedArrhenius(x, y)
```

You can obtain the values of Arrhenius parameters A [1/min^n] and Ea [J/mol]
  ```Python
arrh.arrhenius_params()
```
  ```
Out
{'Ea [J/mol]': 18300.09180621676, 
'A [same as k]': 0.008242440438229868}
```
You can also obtain the plot and assess the fitting.
  ```Python
arrh.plot_arrhenius_fit()
```
![arrhenius](https://user-images.githubusercontent.com/91277572/208984298-0a54b650-b829-4062-8aaf-3f9bedb324e5.png)

  ```Python
arrh.assess_fit()
```
  ```
Out
{'Arrhenius R2': 0.9967098965071525}
```

#References

[Georgiadis, A. G., Charisiou, N. D., Gaber, S., Polychronopoulou, K., Yentekakis, I. V., & Goula, M. A. (2021). Adsorption of hydrogen sulfide at low temperatures using an industrial molecular sieve: an experimental and theoretical study. Acs Omega, 6(23), 14774-14787.](https://pubs.acs.org/doi/full/10.1021/acsomega.0c06157)

[Siakavelas, G. I., Georgiadis, A. G., Charisiou, N. D., Yentekakis, I. V., & Goula, M. A. (2021). Cost‐Effective Adsorption of Oxidative Coupling‐Derived Ethylene Using a Molecular Sieve. Chemical Engineering & Technology, 44(11), 2041-2048.](https://onlinelibrary.wiley.com/doi/abs/10.1002/ceat.202100147)

