from adsorption.models import (Isotherms, Kinetics, ModifiedArrhenius, AdsorptionDynamics, AdsorptionEnthalpy, IsostericHeat)
import matplotlib.pyplot as plt

kinetic = Kinetics()

kinetic.set_inlet()

kinetic.plot_pfo_fit()

plt.show()