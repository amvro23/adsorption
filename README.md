# adsorption
A python package with the purpose of assessing the most reliable isotherm models (i.e., Langmuir, Freundlich, Temkin, Toth, Sips, DR), kinetic models (i.e., PFO, PSO, Weber-Morris, Avrami, Bangham, Elovich) and Arrhenius parameters (i.e., Ea, A), adsorption dynamic models (i.e., Thomas, Yoon-Nelson, Adams-Bohart), and adsorption enthalpy and entropy.

[Install](#Install) / [Usage](#Usage) /  [Isotherms](#Isotherms) / [Kinetics](#Kinetics) / [Arrhenius](#Arrhenius) / [AdsorptionDynamics](#AdsorptionDynamics) / [AdsorptionEnthalpy](#AdsorptionEnthalpy) / [References](#References) / [Contact](#Contact)

# Install
First, make sure you have a Python 3 environment installed.
To install from github:
```Python
pip install -e git+https://github.com/amvro23/adsorption/#egg=adsorption
```
Note: It might be useful to write "git+https://github.com/amvro23/adsorption/#egg=adsorption" if installing directly from a Python interpreter as # can be interpreted as a comment.

# Usage

```Python
from adsorption.models import (Isotherms, Kinetics, ModifiedArrhenius, AdsorptionDynamics, AdsorptionEnthalpy)
import numpy as np
import pandas as pd
```

# Isotherms

Create an instance. Default optional values are P = 1 atm, Mr = 44.01 g/mol for CO2 adsorption, T = 298.15 K, R = 8.205e-5 atm.m3/mol/K).
```Python
isotherm = Isotherms()
```
Set inlet values for adsorption isotherm equations. Optional parameter x represents the dimensionless equilibrium concentration [%] and optional parameter y represents the equilibrium capacity [mg/g].
```Python
isotherm.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following).
```Python
x = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]) 
y = np.array([52.64, 72.59, 99.67, 141.79, 182.48, 203.68, 203.56, 204.33, 204.90])
isotherm.set_inlet(x, y)
```
You can obtain either a single isotherm model plot (e.g., Langmuir model, Freundlich model),
```Python
isotherm.plot_langmuir_fit()
```
![langmuir](https://user-images.githubusercontent.com/91277572/212330579-c1fedef5-444a-4712-98a7-b8ed61e2dbf7.png)

```Python
isotherm.plot_freundlich_fit()
```
![freundlich](https://user-images.githubusercontent.com/91277572/212330713-9117b5e4-7ba2-462e-bb01-6a09edc09ad0.png)

or you can obtain all isotherm models.
```Python
isotherm.plot_all_models()
```
![all_isotherms](https://user-images.githubusercontent.com/91277572/212330912-1888f094-28f6-44f9-bb0c-146cc679b23f.png)

You can assess the fit of isotherm models.
```Python
isotherm.assess_fit()
```
```
Out
{'Langmuir R2': 0.9817163923417866,
 'Freundlich R2': 0.7463363496394846,
 'Temkin R2': 0.8590282631439854,
 'Toth R2': 0.9786970833618412,
 'Sips R2': 0.9794054428262444,
 'DR R2': 0.7707376307818574}
```
 The package can directly return the best isotherm model.
```Python
isotherm.best_fit()
```
```
Out
The best model is that of Langmuir R2 = 0.9817163923417866
```
You can also obtain a dataframe with all the calculated parameters of isotherm equations with their units,
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
or you can have access in isotherm parameters individually in order to create your own plots.
```Python
isotherm.dubinin_radushkevich_params()
```
```
Out
{'E_DR [J/mol]': 147.67076438840127, 
'qmax_DR [mg/g]': 192.6633619518013}
```
```Python
isotherm.sips_curve(isotherm.x)
```
```
Out
array([ 0.04084583,  0.09719205,  0.18721156,  0.36046663,  0.85595811,
        1.64334054,  3.14432104,  7.33485598, 13.69841098])
```
You can get an excel file with the dataframe with all the calculated parameters.
```Python
isotherm.to_excel("Isotherms")
```

# Kinetics

Create an instance.
```Python
kinetic = Kinetics()
```
Set inlet values for adsorption kinetic equations. Optional parameter x represents the adsorption time [min] and optional parameter y represents the cumulative adsorption capacity [mg/g].
```Python
kinetic.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following). 

The csv file 'adsorption_kinetics.csv' can be found in the following link:

https://github.com/amvro23/essentials-ChEng/tree/master/Adsorption%20Processes
```Python
df_kin = pd.read_csv('adsorption_kinetics.csv')
x = df_kin.loc[:, 'minutes'].values
y = df_kin.loc[:, 'qt'].values
kinetic.set_inlet(x, y)
```
You can obtain either a single kinetic model plot (e.g., Bangham model),
```Python
kinetic.plot_bangham_fit()
```
![bangham](https://user-images.githubusercontent.com/91277572/212331223-b793d59e-ca09-43f4-9073-8276a9046582.png)

or you can obtain all kinetic models.
```Python
kinetic.plot_all_models()
```
![kinetics_all](https://user-images.githubusercontent.com/91277572/212331580-23b4c03d-374d-4f8a-9013-e8ecb94758e7.png)

You can assess the fit of kinetic models.
```Python
kinetic.assess_fit()
```
```
Out
{'PFO R2': 0.9835043733020796,
 'PSO R2': 0.9762472082928675,
 'WEBER-MORRIS R2': 0.9280676006898012,
 'AVRAMI R2': 0.9964111275988798,
 'BANGHAM R2': 0.9980776118158825,
 'ELOVICH R2': 0.874909009634501}
 ```
 The package can directly return the best kinetic model.
 ```Python
kinetic.best_fit()
```
 ```
Out
 The best model is that of BANGHAM R2 = 0.9980776118158825
 ```
You can also obtain a dataframe with all the calculated parameters of kinetic equations with their units,
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
10      qmax_bangham [mg/g]       63.29831524175902
11                n_bangham      1.4416628511124052
12     a_elovich [mg/g/min]       4.602398573502221
13         b_elovich [g/mg]     0.05477589170964753
 ```
or you can have access in kinetic parameters individually in order to create your own plots.
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
kinetic.pso_curve(kinetic.x)
```
  ```
  Out
array([ 0.,  0.08073146,  0.16136812, ..., 67.24288195,
       67.25473634, 67.266584  ])
```
You can get an excel file with the dataframe with all the calculated parameters.
```Python
kinetic.to_excel("Kinetics")
```

# Arrhenius

After finding the best kinetic model, the equation is applied to the data obtained at 3 to 4 different temperatures (at least) in order to create the following instance. In this example Bangham model was the best to describe the process of CO2 adsorption.
  ```Python
arrh = ModifiedArrhenius()
```
Set inlet values for modified Arrhenius equation. Optional parameter x represents the temperatures [K] in which adsorption tests were performed and y represents the k value of the best kinetic model, which also carries the units (e.g., if Bangham model is more suitable k carries the units of k in bangham model) at the corresponding temperatures.
  ```Python
arrh.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following).

```Python
x = np.array([298.15, 308.15, 323.15, 373.15])
y = np.array([0.00478, 0.00583, 0.00728, 0.01956])
arrh.set_inlet(x, y)
```
You can obtain the values of Arrhenius parameters A and Ea. Given that bangham model is more suitable to describe the process, A adopts the units from k_bangham [1/min^n].
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
![arrhenius](https://user-images.githubusercontent.com/91277572/212331978-eeaf215a-9b78-4c47-b7ad-784e1bea20aa.png)

  ```Python
arrh.assess_fit()
```
  ```
Out
{'Arrhenius R2': 0.9949596964116147}
```

# AdsorptionDynamics

Create an instance. Default optional values are C=0.1, Mr=44.01 g/mol, T=298.15 K, P=1 atm, h=2 cm, r=0.45 cm, Q=100 ml/min, W=1 g, U=0.1, R=8.205e-5 atm.m3/mol/K) where C represents the initial concentration of the adsorbed molecule (CO2: 10%).
```Python
ads_dyn = AdsorptionDynamics()
```
Set inlet values for adsorption dynamic equations. Optional parameter x represents the adsorption time [min] and optional parameter y represents the dimensionless concentration Ct/C0.
```Python
ads_dyn.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following).

The csv file ('Co_10%.csv') can be found in the following link:

https://github.com/amvro23/essentials-ChEng/tree/master/Adsorption%20Processes
```Python
df_dyn = pd.read_csv('Co_10%.csv')
x = df_dyn.loc[:, 'x'].values
y = df_dyn.loc[:, 'y'].values
ads_dyn.set_inlet(x, y)
```
You can obtain a single adsorption dynamic model plot (e.g., Yoon-nelson model).
```Python
ads_dyn.plot_yoon_nelson_fit()
```
![yoon_nelson](https://user-images.githubusercontent.com/91277572/212332235-925a670d-e7dd-4cc0-9afa-f7eec3f99ba3.png)

Note that Adams-Bohart fit is used for the description of the initial part of the breakthrough curve ranging from 10-15% (default value is U=0.1).

```Python
ads_dyn.plot_adams_bohart_fit()
```
![adams_bohart](https://user-images.githubusercontent.com/91277572/212332359-5aae8a1e-c02b-4c2a-b0a7-629b792ec550.png)

You can assess the fit of dynamic models.
```Python
ads_dyn.assess_fit()
```
```
Out
{'THOMAS R2': 0.9958327653626917,
 'YOON-NELSON R2': 0.9958327653619077,
 'ADAMS-BOHART R2': 0.9983541738775606}
```
You can also obtain a dataframe with all the calculated parameters of dynamic adsorption equations with their units.
```Python
ads_dyn.all_params()
```
```
Out
                          Parameters               Values
0               k_thomas [ml/mg/min]    2.678369108434965
1                 qmax_thomas [mg/g]   202.16211953564525
2              k_yoon_nelson [1/min]   0.4818450777710938
3              tau_yoon_nelson [min]   11.237299942413042
4         k_adams_bohart [ml/mg/min]    5.913583287817446
5            N0_adams_bohart [mg/ml]    137.6058956768907
```
You can get an excel file with the dataframe with all the calculated parameters.
```Python
ads_dyn.to_excel("AdsDynamics")
```

# AdsorptionEnthalpy

Create an instance. Default optional values are C=0.01, Mr=44.01 g/mol, T=298.15 K, P=1 atm, R=8.205e-5 atm.m3/mol/K.

```Python
ads_H = AdsorptionEnthalpy()
```
Set inlet values for adsorption kinetics. Optional parameter x represents the adsorption temperature [K] and optional parameter y represents the equilibrium adsorption capacity [mg/g].
```Python
ads_H.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following).
```Python
x = np.array([298.15, 308.15, 323.15, 348.15, 373.15])
y = np.array([203.6870035, 162.2365645, 116.2852302, 65.14332759, 34.46486588])
ads_H.set_inlet(x, y)
```
You can obtain the values of Vant Hoff parameters enthalpy and entropy,
```Python
ads_H.vant_hoff_params()
```
```
Out
{'enthalpy [J/mol]': -21754.33170685983,
 'entropy [J/mol/K]': -52.31720971406195,
 'R2': 0.9923449042417911,
 'slope': 2616.5902943059696,
 'intercept': -6.292664146507331}
```

```Python
ads_H.plot_vant_hoff()
```
![enthalpy](https://user-images.githubusercontent.com/91277572/212332811-dbbf7adf-8918-4907-b8cc-485135d77f8c.png)

or your can obtain the numerical values to create your own plot.
```Python
ads_H.lnKd
```
```
Out
array([2.4267528, 2.19922383, 1.86621433, 1.28675816, 0.65010871])
```
```Python
ads_H.vant_hoff_line(ads_H.x)
```
```
Out
array([2.4834227 , 2.19862352, 1.80447432, 1.22303396, 0.71950333])
```

# References

[Georgiadis, A. G., Charisiou, N. D., Gaber, S., Polychronopoulou, K., Yentekakis, I. V., & Goula, M. A. (2021). Adsorption of hydrogen sulfide at low temperatures using an industrial molecular sieve: an experimental and theoretical study. Acs Omega, 6(23), 14774-14787.](https://pubs.acs.org/doi/full/10.1021/acsomega.0c06157)

[Georgiadis, A. G., Charisiou, N. D., & Goula, M. A. (2020). Removal of hydrogen sulfide from various industrial gases: A review of the most promising adsorbing materials. Catalysts, 10(5), 521.](https://www.mdpi.com/2073-4344/10/5/521)

[Georgiadis, A. G., Charisiou, N., Yentekakis, I. V., & Goula, M. A. (2020). Hydrogen sulfide (H2S) removal via MOFs. Materials, 13(16), 3640.](https://www.mdpi.com/1996-1944/13/16/3640)

[Siakavelas, G. I., Georgiadis, A. G., Charisiou, N. D., Yentekakis, I. V., & Goula, M. A. (2021). Cost‐Effective Adsorption of Oxidative Coupling‐Derived Ethylene Using a Molecular Sieve. Chemical Engineering & Technology, 44(11), 2041-2048.](https://onlinelibrary.wiley.com/doi/abs/10.1002/ceat.202100147)

# Contact
amvro23@gmail.com

