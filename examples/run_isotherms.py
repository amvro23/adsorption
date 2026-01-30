from adsorption import Isotherms
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


iso = Isotherms()

x = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]) 
y = np.array([52.64, 72.59, 99.67, 141.79, 182.48, 203.68, 203.56, 204.33, 204.90])

iso.set_inlet(x, y)

iso.plot_langmuir_fit()

iso.plot_all_models()

print("")
df = pd.DataFrame.from_dict(iso.assess_fit(), orient="index", columns=["R²"])
print(df)
print("")

iso.best_fit()

print("")
print(iso.all_params())
print("")

plt.show()

