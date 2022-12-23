import numpy as np
import pandas as pd

# ------------------------------------------------------------------------------------------------
# ISOTHERM DATA
# ------------------------------------------------------------------------------------------------

x_iso = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])   
y_iso = np.array([52.64, 72.59, 99.67, 141.79, 182.48, 203.68, 203.56, 204.33, 204.90])

# ------------------------------------------------------------------------------------------------
# KINETIC DATA
# ------------------------------------------------------------------------------------------------

df_kin = pd.read_csv('adsorption_kinetics.csv')
x_kin = df_kin.loc[:, 'minutes'].values
y_kin = df_kin.loc[:, 'qt'].values

# ------------------------------------------------------------------------------------------------
# ARRHENIUS DATA
# ------------------------------------------------------------------------------------------------

x_arrh = np.array([298.15, 308.15, 323.15, 373.15])
y_arrh = np.array([0.00478, 0.00583, 0.00728, 0.01956])

# ------------------------------------------------------------------------------------------------
# DYNAMIC DATA
# ------------------------------------------------------------------------------------------------

df_dyn = pd.read_csv('Co_10%.csv')
x_dyn = df_dyn.loc[:, 'x'].values
y_dyn = df_dyn.loc[:, 'y'].values

# ------------------------------------------------------------------------------------------------
# ENTHALPY DATA
# ------------------------------------------------------------------------------------------------

x_h = np.array([298.15, 308.15, 323.15, 348.15, 373.15])
y_h = np.array([203.6870035, 162.2365645, 116.2852302, 65.14332759, 34.46486588])