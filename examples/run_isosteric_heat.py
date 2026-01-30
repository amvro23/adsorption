from adsorption import IsostericHeat
from adsorption import x_iheat1, x_iheat2, y_iheat1, y_iheat2
import matplotlib.pyplot as plt
import pandas as pd


iso_heat = IsostericHeat()

iso_heat.set_inlet(x1=x_iheat1, y1=y_iheat1, x2=x_iheat2, y2=y_iheat2)

iso_heat.plot_freundlich_langmuir()
plt.show()

print("")
df = pd.DataFrame.from_dict(iso_heat.assess_fit(), orient="index", columns=["R²"])
print(df)
print("")

print("")
print(iso_heat.all_params())
print("")
