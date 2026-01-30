from adsorption import AdsorptionDynamics
from adsorption import x_dyn, y_dyn
import matplotlib.pyplot as plt
import pandas as pd

ads_dyn = AdsorptionDynamics()

ads_dyn.set_inlet(x_dyn, y_dyn)

ads_dyn.plot_thomas_fit()
ads_dyn.plot_yoon_nelson_fit()
ads_dyn.plot_adams_bohart_fit()

print("")
df = pd.DataFrame.from_dict(ads_dyn.assess_fit(), orient="index", columns=["R²"])
print(df)
print("")

ads_dyn.best_fit()

print("")
print(ads_dyn.all_params())
print("")

plt.show()