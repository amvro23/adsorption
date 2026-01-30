from adsorption import ModifiedArrhenius
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

arrh = ModifiedArrhenius()

x = np.array([298.15, 308.15, 323.15, 373.15])
y = np.array([0.00478, 0.00583, 0.00728, 0.01956])

arrh.set_inlet(x, y)

arrh.plot_arrhenius_fit()

print("")
df = pd.DataFrame.from_dict(arrh.arrhenius_params(), orient="index", columns=["Value"])
print(df)
print("")

plt.show()
