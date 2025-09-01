<img width="300" height="300" alt="logo_adsorption" src="https://github.com/user-attachments/assets/23004156-62d1-4533-85b8-e5b76bba37ec" />


# adsorption
A python package with the purpose of assessing the most reliable isotherm models (i.e., Langmuir, Freundlich, Temkin, Toth, Sips, DR), kinetic models (i.e., PFO, PSO, Weber-Morris, Avrami, Bangham, Elovich) and Arrhenius parameters (i.e., Ea, A), adsorption dynamic models (i.e., Thomas, Yoon-Nelson, Adams-Bohart), adsorption enthalpy and entropy as well as isosteric heat of adsorption. The package also offers the opportunity of calculating the amount of adsorbent needed in a scaled-up adsorption unit.

[Install](#Install) / [Usage](#Usage) /  [Isotherms](#Isotherms) / [Kinetics](#Kinetics) / [Arrhenius](#Arrhenius) / [AdsorptionDynamics](#AdsorptionDynamics) / [AdsorptionEnthalpy](#AdsorptionEnthalpy) / [IsostericHeat](#IsostericHeat) / [ScaleUP](#ScaleUP) / [References](#References) / [Contact](#Contact)

# Install
First, make sure you have a Python 3 environment installed.
To install from github:
```Python
pip install -e git+https://github.com/amvro23/adsorption/#egg=adsorption
```
Note: It might be useful to write "git+https://github.com/amvro23/adsorption/#egg=adsorption" if installing directly from a Python interpreter as # can be interpreted as a comment.

# Usage

```Python
from adsorption.models import (Isotherms, Kinetics, ModifiedArrhenius, 
                               AdsorptionDynamics, AdsorptionEnthalpy, IsostericHeat, Adsorbent_ScaleUp)
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
You can also have access to the predicted y-values of each model (e.g., Langmuir model).
```Python
x_obs = iso.x_obs # experimental data

x_fit = np.linspace(min(x_obs), max(x_obs), 100) # x values for plotting the fit
y_fit = iso.langmuir_curve(x_fit) # y values for plotting the fit
```
```
Out
array([ 37.46368625, 119.6448883 , 149.17682327, 164.37925491,
       173.645366  , 179.88424464, 184.37108119, 187.75303548,
       190.39348244, 192.51219663, 194.24991003, 195.70091549,
       196.93076205, 197.98643204, 198.9024801 , 199.7048875 ,
       200.41356577, 201.04403141, 201.60855509, 202.11696822,
       202.57724036, 202.99590002, 203.37834633, 203.72908308,
       204.05189698, 204.34999486, 204.62611061, 204.88258925,
       205.12145371, 205.34445828, 205.55313178, 205.74881266,
       205.93267772, 206.1057658 , 206.26899745, 206.42319126,
       206.56907762, 206.7073103 , 206.83847619, 206.96310371,
       207.08166984, 207.19460631, 207.30230475, 207.40512133,
       207.50338062, 207.59737905, 207.68738789, 207.77365586,
       207.85641142, 207.93586481, 208.0122098 , 208.08562531,
       208.15627678, 208.22431742, 208.28988932, 208.35312445,
       208.41414553, 208.47306686, 208.52999496, 208.58502931,
       208.63826285, 208.68978255, 208.73966985, 208.78800114,
       208.8348481 , 208.88027807, 208.92435438, 208.96713665,
       209.00868103, 209.0490405 , 209.08826503, 209.12640183,
       209.16349552, 209.19958834, 209.23472023, 209.26892909,
       209.3022508 , 209.33471945, 209.36636737, 209.39722529,
       209.42732244, 209.4566866 , 209.48534422, 209.51332051,
       209.54063948, 209.56732401, 209.59339597, 209.61887621,
       209.64378464, 209.66814029, 209.69196137, 209.71526529,
       209.7380687 , 209.76038757, 209.78223716, 209.80363213,
       209.82458652, 209.8451138 , 209.8652269 , 209.88493822])
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
Set inlet values for adsorption kinetic equations. Optional parameter x represents the adsorption time [min] and optional parameter y represents the accumulative adsorption capacity [mg/g].
```Python
kinetic.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following). 

```Python
from adsorption.ads_data import kinetic_data
from io import StringIO

df = pd.read_csv(StringIO(kinetic_data))

# Convert to numpy arrays
x = df["x"].to_numpy()
y = df["y"].to_numpy()

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
![image](https://user-images.githubusercontent.com/91277572/215773255-cc467731-f111-4d9e-8d4a-644f7054cf38.png)

You can assess the fit of kinetic models.
```Python
kinetic.assess_fit()
```
```
Out
{'PFO R2': 0.9835043733020796,
 'PSO R2': 0.9762472082928675,
 'WEBER-MORRIS R2': 0.9280676006898012,
 'AVRAMI R2': 0.9980776118158866,
 'BANGHAM R2': 0.9980776118158825,
 'ELOVICH R2': 0.874909009634501}
 ```
 The package can directly return the best kinetic model.
 ```Python
kinetic.best_fit()
```
 ```
Out
 The best model is that of AVRAMI R2 = 0.9980776118158866
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
6          k_avrami [1/min]    0.024543231342313475
7        qmax_avrami [mg/g]      63.298324486155906
8                  n_avrami      1.4416620072834452
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
{'k_avrami [1/min]': 0.024543231342313475,
 'qmax_avrami [mg/g]': 63.298324486155906,
 'n_avrami': 1.4416620072834452}
```
You can also have access to the predicted values and plot then yourself.
  ```Python
from adsorption.ads_data import kinetic_data
from io import StringIO

df = pd.read_csv(StringIO(kinetic_data))

# Convert to numpy arrays
x = df["x"].to_numpy()
y = df["y"].to_numpy()

kinetic.set_inlet(x, y)
x_fit = np.linspace(min(x), max(x), 100) # x values for plotting the fit
y_fit = kinetic.bangham_curve(x_fit) # y values for plotting the fit

plt.plot(x, y, '.', label='data')
plt.plot(x_fit, y_fit, '--', label='fit')
plt.legend()
```
```
Out
array([ 0.        ,  0.38867841,  1.05022652,  1.87185915,  2.8123119 ,
        3.84648434,  4.95647872,  6.12842533,  7.35102169,  8.61475132,
        9.91142188, 11.23387104, 12.57576851, 13.93147625, 15.29594606,
       16.66464188, 18.03347921, 19.39877653, 20.75721554, 22.1058078 ,
       23.44186635, 24.76298107, 26.06699703, 27.35199516, 28.61627488,
       29.85833828, 31.07687564, 32.27075199, 33.43899477, 34.58078216,
       35.6954323 , 36.78239308, 37.84123255, 38.87162993, 39.87336702,
       40.84632021, 41.79045282, 42.70580794, 43.59250158, 44.45071628,
       45.28069503, 46.0827355 , 46.85718468, 47.60443374, 48.32491327,
       49.01908874, 49.68745628, 50.33053872, 50.94888186, 51.543051  ,
       52.11362771, 52.66120682, 53.1863936 , 53.6898012 , 54.17204822,
       54.63375652, 55.07554918, 55.49804861, 55.90187489, 56.28764419,
       56.65596737, 57.00744868, 57.34268469, 57.66226318, 57.96676229,
       58.25674972, 58.532782  , 58.79540392, 59.04514797, 59.28253397,
       59.50806865, 59.72224539, 59.925544  , 60.11843058, 60.30135734,
       60.47476265, 60.63907095, 60.79469282, 60.94202505, 61.08145077,
       61.21333952, 61.33804754, 61.45591783, 61.56728047, 61.67245281,
       61.77173974, 61.86543393, 61.95381616, 62.03715557, 62.11570996,
       62.18972616, 62.25944024, 62.32507791, 62.38685481, 62.44497683,
       62.49964041, 62.55103289, 62.59933283, 62.64471028, 62.68732713])
```
You can get an excel file with the dataframe with all the calculated parameters.
```Python
kinetic.to_excel("Kinetics")
```

# Arrhenius

After finding the best kinetic model, the equation is applied to the data obtained at 3 to 4 different temperatures (at least) in order to create the following instance. In this example Bangham model is regarded as the best to describe the process of CO2 adsorption.
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
Set inlet values for adsorption enthalpy. Optional parameter x represents the adsorption temperature [K] and optional parameter y represents the equilibrium adsorption capacity [mg/g].
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
Tha package returns directly the vant Hoff plot,
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

# IsostericHeat

Create an instance. Default optional values are T1=273.15 K, T2=293.15 K.
```Python
iso_heat = IsostericHeat()
```
Set inlet values for adsorption isosteric heat. Optional parameter x1 represents the absolute pressure at kPa for adsorption isotherm at T1 and y1 represents the equilibrium adsorption capacity [mmol/g] for adsorption isotherm at T1. The same stands for x2 and y2 but for adsorption isotherm at T2.
```Python
iso_heat.set_inlet()
```
Adjust the values of x and y parameters according to your experimental results (the default values are the following).
```Python
df_iheat1 = pd.read_csv('iso_heat1.csv')
x_iheat1 = df_iheat1.loc[:, 'x'].values
y_iheat1 = df_iheat1.loc[:, 'y'].values

df_iheat2 = pd.read_csv('iso_heat2.csv')
x_iheat2 = df_iheat2.loc[:, 'x'].values
y_iheat2 = df_iheat2.loc[:, 'y'].values

iso_heat.set_inlet(x1=x_iheat1, y1=y_iheat1, x2=x_iheat2, y2=y_iheat2)
```
Tha package returns instantly the Freundlich-Langmuir plot,

```Python
iso_heat.plot_freundlich_langmuir()
```
![image](https://user-images.githubusercontent.com/91277572/215267886-13eb004e-bf9a-4c71-b9d7-2ee4a42b06a4.png)

You can assess the fit of Freundlich-Langmuir model applied at both T1 and T2,
```Python
iso_heat.assess_fit()
```
```
Out
{'Freundlich-Langmuir R2 for T1': 0.9998116123272973,
 'Freundlich-Langmuir R2 for T2': 0.9999895124385948}
```
while you can also get access to the calculated parameters of the freundlich-langmuir model at both T1 and T2.
```Python
iso_heat.all_params()
```
```
Out
                Parameters                Values
0         a at T1 [mmol/g]     5.931481472276362
1        b at T1 [1/kPa^c]   0.05028631551093073
2  c at T1 [dimensionless]    1.0559760862112229
3         a at T2 [mmol/g]     5.892741001110274
4        b at T2 [1/kPa^c]  0.016824160548525938
5  c at T2 [dimensionless]    1.0951835091030404
```

You can also get the dataframe of the desired parameters to create your own plots if necessary,
```Python
iso_heat.get_dataframe()
```
and save them into an excel file.
```Python
iso_heat.to_excel("IsoHeat")
```
The package directly returns the plot of isosteric heat (calculated via the Clausius-Clapeyron equation) vs adsorbed quantity.
```Python
iso_heat.plot_isoHeat_vs_mmol()
```
![image](https://user-images.githubusercontent.com/91277572/215267368-5e9a62c5-f964-48b2-924d-b086300db94e.png)

# ScaleUp

The package offers the opportunity of calculating the amount of adsorbent needed in a scaled-up adsorption unit.

Create an instance. Default optional values are Mr=34.1 as H2S is the targe molecule, and molar fraction is y=0.001.
```Python
scale_up = Adsorbent_ScaleUp()
```

Set inlet values for the experimental and pilot unit. Default values are given herein. Adsorption capacity in mg/g for y=0.001, adsorption time in min for y=0.001, and flow rate in mL/min.
```Python
scale_up.exp_unit(ads_capacity=125, exp_ads_time=500)
scale_up.pilot_unit(total_flow_rate=2500, pilot_ads_time=2000)
```
The amount of the adsorbent needed for the pilot unit is automaticaly printed in grams, however you can also have access to it as an attribute.
```Python
scale_up.quantity_ads_in_pilot_unit
```
```
Out
60.89285714285715
```

# References

[Georgiadis, A. G., Charisiou, N. D., Gaber, S., Polychronopoulou, K., Yentekakis, I. V., & Goula, M. A. (2021). Adsorption of hydrogen sulfide at low temperatures using an industrial molecular sieve: an experimental and theoretical study. Acs Omega, 6(23), 14774-14787.](https://pubs.acs.org/doi/full/10.1021/acsomega.0c06157)

[Georgiadis, A. G., Charisiou, N. D., & Goula, M. A. (2020). Removal of hydrogen sulfide from various industrial gases: A review of the most promising adsorbing materials. Catalysts, 10(5), 521.](https://www.mdpi.com/2073-4344/10/5/521)

[Georgiadis, A. G., Charisiou, N., Yentekakis, I. V., & Goula, M. A. (2020). Hydrogen sulfide (H2S) removal via MOFs. Materials, 13(16), 3640.](https://www.mdpi.com/1996-1944/13/16/3640)

[Siakavelas, G. I., Georgiadis, A. G., Charisiou, N. D., Yentekakis, I. V., & Goula, M. A. (2021). Cost‐Effective Adsorption of Oxidative Coupling‐Derived Ethylene Using a Molecular Sieve. Chemical Engineering & Technology, 44(11), 2041-2048.](https://onlinelibrary.wiley.com/doi/abs/10.1002/ceat.202100147)

[Nuhnena, A., & Christoph Janiak, C. (2020), A practical guide to calculate the isosteric heat/enthalpy of adsorption via adsorption isotherms in metal-organic frameworks, MOFs](https://pubs.rsc.org/en/content/articlelanding/2020/dt/d0dt01784a)

# Contact
amvro23@gmail.com

