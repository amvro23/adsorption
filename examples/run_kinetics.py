from adsorption import Kinetics
from adsorption import kinetic_data
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO


kin = Kinetics()

df = pd.read_csv(StringIO(kinetic_data))

# Convert to numpy arrays
x = df["x"].to_numpy()
y = df["y"].to_numpy()

kin.set_inlet(x, y)

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
