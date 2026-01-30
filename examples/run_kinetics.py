from adsorption import Kinetics
from adsorption import x_kin, y_kin
import matplotlib.pyplot as plt
import pandas as pd


kin = Kinetics()


kin.set_inlet(x_kin, y_kin)

kin.plot_bangham_fit()
kin.plot_all_models()

print("")
df = pd.DataFrame.from_dict(kin.assess_fit(), orient="index", columns=["R²"])
print(df)
print("")

kin.best_fit()

print("")
print(kin.all_params())
print("")

plt.show()
