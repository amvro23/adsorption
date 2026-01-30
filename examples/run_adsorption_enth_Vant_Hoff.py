from adsorption import AdsorptionEnthalpy
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ads_H = AdsorptionEnthalpy()

x = np.array([298.15, 308.15, 323.15, 348.15, 373.15])
y = np.array([203.6870035, 162.2365645, 116.2852302, 65.14332759, 34.46486588])

ads_H.set_inlet(x, y)

ads_H.plot_vant_hoff()

print("")
df = pd.DataFrame.from_dict(ads_H.vant_hoff_params(), orient="index", columns=["Value"])
print(df)
print("")

plt.show()