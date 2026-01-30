import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from .ads_data import (x_iso, y_iso, x_kin, y_kin, x_arrh, y_arrh, 
                             x_dyn, y_dyn, x_h, y_h, x_iheat1, y_iheat1, x_iheat2, y_iheat2)


class Isotherms(object):

    def __init__(self, P=1, Mr=44.01, T=298.15, R=8.205e-5):
        """Class for evaluating 6 different adsorption isotherms."""
        self.P = P
        self.Mr = Mr
        self.T = T
        self.R = R

    def set_inlet(self, x=x_iso, y=y_iso):
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y
        self.factor = (self.P * self.Mr) / self.R / self.T
        self.x_obs = self.factor * self.x
        self.xfit = np.linspace(min(self.x_obs), max(self.x_obs), 50)

    # -------------------------
    # Helpers (private)
    # -------------------------
    def _fit_and_pack(self, model_func, p0, bounds, attr_names, keys):
        FitParams, _ = curve_fit(
            model_func,
            self.x_obs,
            self.y,
            np.array(p0),
            bounds=bounds,
        )
        for name, val in zip(attr_names, FitParams):
            setattr(self, name, val)
        return {k: v for k, v in zip(keys, FitParams)}

    def _curve_from_params(self, model_func, x, params_dict, keys):
        p = [params_dict[k] for k in keys]
        return model_func(x, *p)

    def _plot_fit(self, title, y_pred, x_label="$C_e$ $[mg/L]$", y_label="$q_e$ $[mg/g]$"):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(self.x_obs, self.y, "ko", mfc="none", label="Observed")
        ax.plot(self.xfit, y_pred, "k--", label="Predicted")
        ax.set_xlabel(x_label, fontsize=10, fontweight="bold")
        ax.set_ylabel(y_label, fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title(title)
        ax.grid(ls=":")
        fig.tight_layout()
        return fig, ax

    @staticmethod
    def _adj_r2(y_true, y_pred, n_params):
        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        n = len(y_true)
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        r2 = 1 - ss_res / ss_tot
        return r2 * ((n - 1) / (n - n_params))

    # -------------------------
    # Models (same names)
    # -------------------------
    def langmuir(self, x, k, q):
        x = np.array(x)
        return k * q * x / (1 + k * x)

    def langmuir_params(self):
        keys = ["K_Langmuir [L/mg]", "qmax_Langmuir [mg/g]"]
        return self._fit_and_pack(
            model_func=self.langmuir,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "q"],
            keys=keys,
        )

    def langmuir_curve(self, x):
        keys = ["K_Langmuir [L/mg]", "qmax_Langmuir [mg/g]"]
        return self._curve_from_params(self.langmuir, x, self.langmuir_params(), keys)

    def plot_langmuir_fit(self):
        return self._plot_fit("Langmuir fit", self.langmuir_curve(self.xfit))

    def freundlich(self, x, k, n):
        x = np.array(x)
        return k * x ** (1 / n)

    def freundlich_params(self):
        keys = ["K_Freundlich [L^(1/n_Freundlich)·mg^(1-1/n_Freundlich)·g^-1]", "n_Freundlich"]
        return self._fit_and_pack(
            model_func=self.freundlich,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "n"],
            keys=keys,
        )

    def freundlich_curve(self, x):
        keys = ["K_Freundlich [L^(1/n_Freundlich)·mg^(1-1/n_Freundlich)·g^-1]", "n_Freundlich"]
        return self._curve_from_params(self.freundlich, x, self.freundlich_params(), keys)

    def plot_freundlich_fit(self):
        return self._plot_fit("Freundlich fit", self.freundlich_curve(self.xfit))

    def temkin(self, x, a, b):
        x = np.array(x)
        return (8.314 * 298.15 / b) * (np.log(a * x))

    def temkin_params(self):
        keys = ["A_Temkin [L/mg]", "B_Temkin [J/mol]"]
        return self._fit_and_pack(
            model_func=self.temkin,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["a", "b"],
            keys=keys,
        )

    def temkin_curve(self, x):
        keys = ["A_Temkin [L/mg]", "B_Temkin [J/mol]"]
        return self._curve_from_params(self.temkin, x, self.temkin_params(), keys)

    def plot_temkin_fit(self):
        return self._plot_fit("Temkin fit", self.temkin_curve(self.xfit))

    def toth(self, x, k, n, q):
        x = np.array(x)
        return q * x / (k + x ** n) ** (1 / n)

    def toth_params(self):
        keys = ["K_Toth [mg^z·L^-n_toth]", "n_Toth", "qmax_Toth [mg/g]"]
        return self._fit_and_pack(
            model_func=self.toth,
            p0=[1, 1, 1],
            bounds=(0, [np.inf, np.inf, np.inf]),
            attr_names=["k", "n", "q"],
            keys=keys,
        )

    def toth_curve(self, x):
        keys = ["K_Toth [mg^z·L^-n_toth]", "n_Toth", "qmax_Toth [mg/g]"]
        return self._curve_from_params(self.toth, x, self.toth_params(), keys)

    def plot_toth_fit(self):
        return self._plot_fit("Toth fit", self.toth_curve(self.xfit))

    def sips(self, x, k, n, q):
        x = np.array(x)
        return (q * k * x ** (1 / n)) / (1 + k * x ** (1 / n))

    def sips_params(self):
        keys = ["K_Sips [L^n_Sips·mg^-n_Sips]", "n_Sips", "qmax_Sips [mg/g]"]
        return self._fit_and_pack(
            model_func=self.sips,
            p0=[1, 1, 1],
            bounds=(0, [np.inf, np.inf, np.inf]),
            attr_names=["k", "n", "q"],
            keys=keys,
        )

    def sips_curve(self, x):
        keys = ["K_Sips [L^n_Sips·mg^-n_Sips]", "n_Sips", "qmax_Sips [mg/g]"]
        return self._curve_from_params(self.sips, x, self.sips_params(), keys)

    def plot_sips_fit(self):
        return self._plot_fit("Sips fit", self.sips_curve(self.xfit))

    def dubinin_radushkevich(self, x, E, q):
        x = np.array(x)
        return q * np.exp((8.314 * 298.15 * np.log(1 + 1 / x)) ** 2 / (-2 * E ** 2))

    def dubinin_radushkevich_params(self):
        keys = ["E_DR [J/mol]", "qmax_DR [mg/g]"]
        return self._fit_and_pack(
            model_func=self.dubinin_radushkevich,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["E", "q"],
            keys=keys,
        )

    def dr_curve(self, x):
        keys = ["E_DR [J/mol]", "qmax_DR [mg/g]"]
        return self._curve_from_params(self.dubinin_radushkevich, x, self.dubinin_radushkevich_params(), keys)

    def plot_dr_fit(self):
        return self._plot_fit("DR fit", self.dr_curve(self.xfit))

    def plot_all_models(self):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(self.x_obs, self.y, "ko", mfc="none", label="Observed")
        ax.plot(self.xfit, self.langmuir_curve(self.xfit), "r--", label="Langmuir")
        ax.plot(self.xfit, self.freundlich_curve(self.xfit), "b--", label="Freundlich")
        ax.plot(self.xfit, self.temkin_curve(self.xfit), "g--", label="Temkin")
        ax.plot(self.xfit, self.toth_curve(self.xfit), "y--", label="Toth")
        ax.plot(self.xfit, self.sips_curve(self.xfit), "m--", label="Sips")
        ax.plot(self.xfit, self.dr_curve(self.xfit), "c--", label="DR")
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize=10, fontweight="bold")
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title("All models")
        ax.grid(ls=":")
        fig.tight_layout()
        return fig, ax

    def assess_fit(self):
        y_obs = self.y

        y_langmuir = self.langmuir_curve(self.x_obs)
        y_freundlich = self.freundlich_curve(self.x_obs)
        y_temkin = self.temkin_curve(self.x_obs)
        y_toth = self.toth_curve(self.x_obs)
        y_sips = self.sips_curve(self.x_obs)
        y_dr = self.dr_curve(self.x_obs)

        return {
            "Langmuir R2": self._adj_r2(y_obs, y_langmuir, n_params=len(self.langmuir_params())),
            "Freundlich R2": self._adj_r2(y_obs, y_freundlich, n_params=len(self.freundlich_params())),
            "Temkin R2": self._adj_r2(y_obs, y_temkin, n_params=len(self.temkin_params())),
            "Toth R2": self._adj_r2(y_obs, y_toth, n_params=len(self.toth_params())),
            "Sips R2": self._adj_r2(y_obs, y_sips, n_params=len(self.sips_params())),
            "DR R2": self._adj_r2(y_obs, y_dr, n_params=len(self.dubinin_radushkevich_params())),
        }

    def best_fit(self):
        model = max(self.assess_fit(), key=self.assess_fit().get)
        value = self.assess_fit().get(model)
        return print("The best model is that of", model, "=", value)

    def all_params(self):
        def get_params(*args):
            return np.vstack(args)

        all_p = get_params(
            list(self.langmuir_params().items()),
            list(self.freundlich_params().items()),
            list(self.temkin_params().items()),
            list(self.toth_params().items()),
            list(self.sips_params().items()),
            list(self.dubinin_radushkevich_params().items()),
        )

        df = pd.DataFrame(all_p, columns=["Parameters", "Values"])
        return df

    def to_excel(self, filename, **options):
        path = filename + ".xlsx"
        with pd.ExcelWriter(path) as writer:
            self.all_params().to_excel(writer, sheet_name="Isotherms")



class Kinetics(object):

    def __init__(self):
        """Class for evaluating 6 different adsorption kinetic models."""

    def set_inlet(self, x=x_kin, y=y_kin):
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y

    # -------------------------
    # Helpers (private)
    # -------------------------
    def _fit_and_pack(self, model_func, xdata, ydata, p0, bounds, attr_names, keys):
        FitParams, _ = curve_fit(model_func, xdata, ydata, np.array(p0), bounds=bounds)

        for name, val in zip(attr_names, FitParams):
            setattr(self, name, val)

        return {k: v for k, v in zip(keys, FitParams)}

    def _curve_from_params(self, model_func, x, params_dict, keys):
        p = [params_dict[k] for k in keys]
        return model_func(x, *p)

    def _plot_fit(self, title, y_pred, x_label="$Time$ $[min]$", y_label="$q_t$ $[mg/g]$", ylim_pad=20):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(self.x, self.y, "k:", mfc="none", label="Observed")
        ax.plot(self.x, y_pred, "r--", label="Predicted")
        ax.set_xlabel(x_label, fontsize=10, fontweight="bold")
        ax.set_ylabel(y_label, fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title(title)
        ax.grid(ls=":")
        if ylim_pad is not None:
            ax.set_ylim(0, np.max(self.y) + ylim_pad)
        fig.tight_layout()
        return fig, ax

    @staticmethod
    def _adj_r2(y_true, y_pred, n_params):
        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        n = len(y_true)
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        r2 = 1 - ss_res / ss_tot
        return r2 * ((n - 1) / (n - n_params))

    # -------------------------
    # Models (unchanged names)
    # -------------------------
    def pfo(self, x, k, q):
        x = np.array(x)
        return q * (1 - np.exp(-k * x))

    def pfo_params(self):
        keys = ["k_pfo [1/min]", "qmax_PFO [mg/g]"]
        return self._fit_and_pack(
            model_func=self.pfo,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "q"],
            keys=keys,
        )

    def pfo_curve(self, x):
        keys = ["k_pfo [1/min]", "qmax_PFO [mg/g]"]
        return self._curve_from_params(self.pfo, x, self.pfo_params(), keys)

    def plot_pfo_fit(self):
        return self._plot_fit("PFO fit", self.pfo_curve(self.x))

    def pso(self, x, k, q):
        x = np.array(x)
        return k * q ** 2 * x / (1 + k * q * x)

    def pso_params(self):
        keys = ["k_pso [g mg^-1 min^-1]", "qmax_PSO [mg/g]"]
        return self._fit_and_pack(
            model_func=self.pso,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "q"],
            keys=keys,
        )

    def pso_curve(self, x):
        keys = ["k_pso [g mg^-1 min^-1]", "qmax_PSO [mg/g]"]
        return self._curve_from_params(self.pso, x, self.pso_params(), keys)

    def plot_pso_fit(self):
        return self._plot_fit("PSO fit", self.pso_curve(self.x))

    def weber_morris(self, x, k, c):
        x = np.array(x)
        return k * x ** 0.5 + c

    def weber_morris_params(self):
        keys = ["k_wm [mg g^-1 min^-0.5]", "C"]
        return self._fit_and_pack(
            model_func=self.weber_morris,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "c"],
            keys=keys,
        )

    def weber_morris_curve(self, x):
        keys = ["k_wm [mg g^-1 min^-0.5]", "C"]
        return self._curve_from_params(self.weber_morris, x, self.weber_morris_params(), keys)

    def plot_weber_morris_fit(self):
        return self._plot_fit("Weber-Morris fit", self.weber_morris_curve(self.x))

    def avrami(self, x, k, q, n):
        x = np.array(x)
        return q * (1 - np.exp(-(k * x) ** n))

    def avrami_params(self):
        keys = ["k_avrami [1/min]", "qmax_avrami [mg/g]", "n_avrami"]
        return self._fit_and_pack(
            model_func=self.avrami,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1, 1],
            bounds=(0, [np.inf, np.inf, np.inf]),
            attr_names=["k", "q", "n"],
            keys=keys,
        )

    def avrami_curve(self, x):
        keys = ["k_avrami [1/min]", "qmax_avrami [mg/g]", "n_avrami"]
        return self._curve_from_params(self.avrami, x, self.avrami_params(), keys)

    def plot_avrami_fit(self):
        return self._plot_fit("Avrami fit", self.avrami_curve(self.x))

    def bangham(self, x, k, q, n):
        x = np.array(x)
        return q * (1 - np.exp(-k * x ** n))

    def bangham_params(self):
        keys = ["k_bangham [1/min^n]", "qmax_bangham [mg/g]", "n_bangham"]
        return self._fit_and_pack(
            model_func=self.bangham,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1, 1],
            bounds=(0, [np.inf, np.inf, np.inf]),
            attr_names=["k", "q", "n"],
            keys=keys,
        )

    def bangham_curve(self, x):
        keys = ["k_bangham [1/min^n]", "qmax_bangham [mg/g]", "n_bangham"]
        return self._curve_from_params(self.bangham, x, self.bangham_params(), keys)

    def plot_bangham_fit(self):
        return self._plot_fit("Bangham fit", self.bangham_curve(self.x))

    def elovich(self, x, a, b):
        x = np.array(x)
        return (1 / b) * np.log(a * b * x)

    def elovich_params(self):
        keys = ["a_elovich [mg/g/min]", "b_elovich [g/mg]"]
        # your original fit uses x[1:], y[1:] to avoid log(0)
        return self._fit_and_pack(
            model_func=self.elovich,
            xdata=self.x[1:],
            ydata=self.y[1:],
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["a", "b"],
            keys=keys,
        )

    def elovich_curve(self, x):
        keys = ["a_elovich [mg/g/min]", "b_elovich [g/mg]"]
        return self._curve_from_params(self.elovich, x, self.elovich_params(), keys)

    def plot_elovich_fit(self):
        # same y-lim behavior as your original
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(self.x, self.y, "k:", mfc="none", label="Observed")
        ax.plot(self.x, self.elovich_curve(self.x), "r--", label="Predicted")
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight="bold")
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title("Elovich fit")
        ax.set_ylim(0, np.max(self.y) + 20)
        ax.grid(ls=":")
        fig.tight_layout()
        return fig, ax

    def plot_all_models(self):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(self.x, self.y, "k:", mfc="none", label="Observed")
        ax.plot(self.x, self.pfo_curve(self.x), "r--", label="PFO")
        ax.plot(self.x, self.pso_curve(self.x), "b--", label="PSO")
        ax.plot(self.x, self.weber_morris_curve(self.x), "g--", label="Weber-Morris")
        ax.plot(self.x, self.avrami_curve(self.x), "y--", label="Avrami")
        ax.plot(self.x, self.bangham_curve(self.x), "m--", label="Bangham")
        ax.plot(self.x, self.elovich_curve(self.x), "c--", label="Elovich")
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight="bold")
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title("All models")
        ax.set_ylim(0, np.max(self.y) + 20)
        ax.grid(ls=":")
        fig.tight_layout()
        return fig, ax

    def assess_fit(self):
        y_observed = self.y

        y_pfo = self.pfo_curve(self.x)
        y_pso = self.pso_curve(self.x)
        y_wm = self.weber_morris_curve(self.x)
        y_avrami = self.avrami_curve(self.x)
        y_bangham = self.bangham_curve(self.x)

        # Elovich was fitted on x[1:], y[1:], so assess similarly
        y_elovich = self.elovich_curve(self.x[1:])

        return {
            "PFO R2": self._adj_r2(y_observed, y_pfo, n_params=len(self.pfo_params())),
            "PSO R2": self._adj_r2(y_observed, y_pso, n_params=len(self.pso_params())),
            "WEBER-MORRIS R2": self._adj_r2(y_observed, y_wm, n_params=len(self.weber_morris_params())),
            "AVRAMI R2": self._adj_r2(y_observed, y_avrami, n_params=len(self.avrami_params())),
            "BANGHAM R2": self._adj_r2(y_observed, y_bangham, n_params=len(self.bangham_params())),
            "ELOVICH R2": self._adj_r2(y_observed[1:], y_elovich, n_params=len(self.elovich_params())),
        }

    def best_fit(self):
        model = max(self.assess_fit(), key=self.assess_fit().get)
        value = self.assess_fit().get(model)
        return print("The best model is that of", model, "=", value)

    def all_params(self):
        def get_params(*args):
            return np.vstack(args)

        all_p = get_params(
            list(self.pfo_params().items()),
            list(self.pso_params().items()),
            list(self.weber_morris_params().items()),
            list(self.avrami_params().items()),
            list(self.bangham_params().items()),
            list(self.elovich_params().items()),
        )

        df = pd.DataFrame(all_p, columns=["Parameters", "Values"])
        return df

    def to_excel(self, filename, **options):
        path = filename + ".xlsx"
        with pd.ExcelWriter(path) as writer:
            self.all_params().to_excel(writer, sheet_name="Kinetics")

    
class ModifiedArrhenius(object):
    
    def __init__(self):
        """Class for calculating Arrhenius parameters. 
        """
        
    def set_inlet(self, x=x_arrh, y=y_arrh):
        """Set inlet parameters for modified Arrhenius equation.
        Parameters
        ----------
        x  : 1d array of floats, optional
             Adsorption temperatures [K]
        y  : 1d array of floats, optional
             Kinetic constants of the best fitted model at different temperatures 
        """
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y
        
    def arrhenius(self, x, Ea, A):
        return A*np.exp(-Ea*(1/x - 1/np.mean(self.x))*(1/8.314))
    
    def arrhenius_params(self):
        FitParams, FitCov = curve_fit(self.arrhenius,
                                      self.x, 
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.Ea = FitParams[0]
        self.A = FitParams[1]
        return {"Ea [J/mol]": self.Ea,
                "A [same as k]": self.A}

    def arrhenius_curve(self, x):
        yfit = self.arrhenius(x, 
                              self.arrhenius_params()["Ea [J/mol]"],
                              self.arrhenius_params()["A [same as k]"])
        return yfit
        
    def plot_arrhenius_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        xfit = np.linspace(min(self.x), max(self.x), 200)
        ax.plot(self.x, self.y, "ko", mfc = "none", label = "Observed")
        ax.plot(xfit, self.arrhenius_curve(xfit), "g--", label = "Arrhenius")
        ax.set_xlabel("$Temperature$ $[K]$", fontsize=10, fontweight="bold")
        ax.set_ylabel("$k$", fontsize=10, fontweight="bold")
        ax.set_title("Modified Arrhenius", fontsize=10)
        ax.grid(linestyle=":")
        fig.tight_layout()
        
    def assess_fit(self):
        y_observed = self.y
        y_arrh = self.arrhenius_curve(self.x)
        n = len(y_observed)
        R_arrhenius = 1-((np.sum((y_observed - y_arrh)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.arrhenius_params().items())))
        return {"Arrhenius R2": R_arrhenius}


class AdsorptionDynamics(object):

    def __init__(self, C=0.1, Mr=44.01, T=298.15,
                 P=1, h=2, r=0.45, Q=100, W=1, U=0.1, R=8.205e-5):
        """Class for calculating 3 different empirical adsorption dynamic models."""
        self.C = C
        self.Mr = Mr
        self.T = T
        self.P = P
        self.h = h
        self.r = r
        self.Q = Q
        self.W = W
        self.U = U
        self.R = R

        self.c0 = self.C * (self.P * self.Mr) / self.R / self.T / 1000
        self.A = np.pi * (self.r ** 2)
        self.v = self.Q / self.A

    def set_inlet(self, x=x_dyn, y=y_dyn):
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y
        self.yy = self.y[self.y <= self.U]
        self.xx = self.x[np.where(self.y <= self.U)[0]]

    # -------------------------
    # Helpers (private)
    # -------------------------
    def _fit_and_pack(self, model_func, xdata, ydata, p0, bounds, attr_names, keys):
        FitParams, _ = curve_fit(
            model_func,
            xdata,
            ydata,
            np.array(p0),
            bounds=bounds,
        )
        for name, val in zip(attr_names, FitParams):
            setattr(self, name, val)
        return {k: v for k, v in zip(keys, FitParams)}

    def _curve_from_params(self, model_func, x, params_dict, keys):
        p = [params_dict[k] for k in keys]
        return model_func(x, *p)

    def _plot_fit(self, title, x_obs, y_obs, x_pred, y_pred,
                  x_label="$Time$ $[min]$", y_label="$C_t/C_0$"):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(x_obs, y_obs, "k:", mfc="none", label="Observed")
        ax.plot(x_pred, y_pred, "r--", label="Predicted")
        ax.set_xlabel(x_label, fontsize=10, fontweight="bold")
        ax.set_ylabel(y_label, fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title(title)
        ax.grid(ls=":")
        fig.tight_layout()
        return fig, ax

    @staticmethod
    def _adj_r2(y_true, y_pred, n_params):
        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        n = len(y_true)
        ss_res = np.sum((y_true - y_pred) ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        r2 = 1 - ss_res / ss_tot
        return r2 * ((n - 1) / (n - n_params))

    # -------------------------
    # Thomas (same names)
    # -------------------------
    def thomas(self, x, k, q):
        x = np.array(x)
        return 1 / (1 + np.exp(k * q * (self.W / self.Q) - k * self.c0 * x))

    def thomas_params(self):
        keys = ["k_thomas [ml/mg/min]", "qmax_thomas [mg/g]"]
        return self._fit_and_pack(
            model_func=self.thomas,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "q"],
            keys=keys,
        )

    def thomas_curve(self, x):
        keys = ["k_thomas [ml/mg/min]", "qmax_thomas [mg/g]"]
        return self._curve_from_params(self.thomas, x, self.thomas_params(), keys)

    def plot_thomas_fit(self):
        return self._plot_fit(
            title="Thomas fit",
            x_obs=self.x, y_obs=self.y,
            x_pred=self.x, y_pred=self.thomas_curve(self.x),
        )

    # -------------------------
    # Yoon-Nelson (same names)
    # -------------------------
    def yoon_nelson(self, x, k, tau):
        x = np.array(x)
        return np.exp(k * (x - tau)) / (1 + np.exp(k * (x - tau)))

    def yoon_nelson_params(self):
        keys = ["k_yoon_nelson [1/min]", "tau_yoon_nelson [min]"]
        return self._fit_and_pack(
            model_func=self.yoon_nelson,
            xdata=self.x,
            ydata=self.y,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "tau"],
            keys=keys,
        )

    def yoon_nelson_curve(self, x):
        keys = ["k_yoon_nelson [1/min]", "tau_yoon_nelson [min]"]
        return self._curve_from_params(self.yoon_nelson, x, self.yoon_nelson_params(), keys)

    def plot_yoon_nelson_fit(self):
        return self._plot_fit(
            title="Yoon-Nelson fit",
            x_obs=self.x, y_obs=self.y,
            x_pred=self.x, y_pred=self.yoon_nelson_curve(self.x),
        )

    # -------------------------
    # Adams-Bohart (same names)
    # -------------------------
    def adams_bohart(self, x, k, N0):
        x = np.array(x)
        return np.exp(k * self.c0 * x - k * N0 * (self.h / self.v))

    def adams_bohart_params(self):
        keys = ["k_adams_bohart [ml/mg/min]", "N0_adams_bohart [mg/ml]"]
        return self._fit_and_pack(
            model_func=self.adams_bohart,
            xdata=self.xx,
            ydata=self.yy,
            p0=[1, 1],
            bounds=(0, [np.inf, np.inf]),
            attr_names=["k", "N0"],
            keys=keys,
        )

    def adams_bohart_curve(self, x):
        keys = ["k_adams_bohart [ml/mg/min]", "N0_adams_bohart [mg/ml]"]
        return self._curve_from_params(self.adams_bohart, x, self.adams_bohart_params(), keys)

    def plot_adams_bohart_fit(self):
        # Observed is full curve (x,y), predicted is only the early-stage region (xx,yy)
        return self._plot_fit(
            title="Adams-Bohart fit",
            x_obs=self.x, y_obs=self.y,
            x_pred=self.xx, y_pred=self.adams_bohart_curve(self.xx),
        )

    # -------------------------
    # Combined
    # -------------------------
    def plot_all_models(self):
        fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        ax.plot(self.x, self.y, "k:", mfc="none", label="Observed")
        ax.plot(self.x, self.thomas_curve(self.x), "r--", label="Thomas")
        ax.plot(self.x, self.yoon_nelson_curve(self.x), "b--", label="Yoon-Nelson")
        ax.plot(self.xx, self.adams_bohart_curve(self.xx), "m--", label="Adams-Bohart")
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight="bold")
        ax.set_ylabel("$C_t/C_0$", fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title("All models")
        ax.grid(ls=":")
        fig.tight_layout()
        return fig, ax

    def assess_fit(self):
        y_obs1 = self.y
        y_obs2 = self.yy

        y_thomas = self.thomas_curve(self.x)
        y_yoon_nelson = self.yoon_nelson_curve(self.x)
        y_adams_bohart = self.adams_bohart_curve(self.xx)

        return {
            "THOMAS R2": self._adj_r2(y_obs1, y_thomas, n_params=len(self.thomas_params())),
            "YOON-NELSON R2": self._adj_r2(y_obs1, y_yoon_nelson, n_params=len(self.yoon_nelson_params())),
            "ADAMS-BOHART R2": self._adj_r2(y_obs2, y_adams_bohart, n_params=len(self.adams_bohart_params())),
        }

    def best_fit(self):
        model = max(self.assess_fit(), key=self.assess_fit().get)
        value = self.assess_fit().get(model)
        return print("The best model is that of", model, "=", value)

    def all_params(self):
        def get_params(*args):
            params = np.vstack(args)
            return params

        all_p = get_params(
            list(self.thomas_params().items()),
            list(self.yoon_nelson_params().items()),
            list(self.adams_bohart_params().items()),
        )

        df = pd.DataFrame(all_p, columns=["Parameters", "Values"])
        return df

    def to_excel(self, filename, **options):
        path = filename + ".xlsx"
        with pd.ExcelWriter(path) as writer:
            self.all_params().to_excel(writer, sheet_name="AdsDynamics")


class AdsorptionEnthalpy(object):
    
    def __init__(self, C=0.01, Mr=44.01, T=298.15, P=1, R=8.205e-5):
        """Class for assessing thermodynamic parameters (i.e., Enthalpy, Entropy).            
        Parameters
        ----------
        C  : float or integer, optional
             Initial concentration of adsorbed molecule [%], by default 0.01 (1% CO2)
        Mr : float, optional
             Molar mass of adsorbed molecule [g/mol], by default 44.01 for CO2
        T  : float or integer, optional
             Adsorption temperature [K], by default 298.15
        P  : float or integer, optional
             Adsorption pressure [atm], by default 1
        R  : float, optional
             Gas constant [atm.m3/mol/K], by default 8.205e-5
        
        Calculated Parameters
        ----------
        c0      : float
                  Converted equilibrium concentration [mg/L]
        lnKd    : float
                  Ratio of ln(qe/Ce)
        x_inv   : inverse predictor
        """   
        self.C = C
        self.Mr = Mr
        self.T = T
        self.P = P
        self.R = R
        self.c0 = self.C*(self.P*self.Mr)/self.R/self.T
        
    def set_inlet(self, x=x_h, y=y_h):
        """Set inlet parameters for Vant Hoff equation.
        Parameters
        ----------
        x  : 1d array of floats, optional
             Adsorption temperature [K]
        y  : 1d array of floats, optional
             Equilibrium adsorption capacity [mg/g]
        
        Calculated Parameters
        ----------
        lnKd    : float
                  Ratio of ln(qe/Ce)
        x_inv   : inverse predictor
        """
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y
        self.x_inv = 1/self.x
        self.lnKd = np.log(self.y/self.c0)
        
    def vant_hoff_params(self):
        slope, intercept, r, p, std_err = stats.linregress(self.x_inv, self.lnKd)
        enthlapy = -slope*8.314
        entropy = intercept*8.314
        return {"enthalpy [J/mol]": enthlapy,
                "entropy [J/mol/K]": entropy,
                "R2": r**2, 
                "slope": slope, 
                "intercept": intercept}
        
    def vant_hoff_line(self, x):
        yfit = list(map(lambda x: self.vant_hoff_params()["slope"]*x 
                        + self.vant_hoff_params()["intercept"], 1/x))
        return np.array(yfit)
    
    def plot_vant_hoff(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_inv, self.lnKd, "ko", mfc = "none", label = "Observed")
        ax.plot(self.x_inv, self.vant_hoff_line(self.x), "r--", label = "Vant Hoff")
        ax.set_xlabel("$1/T$ $[1/K]$", fontsize=10, fontweight="bold")
        ax.set_ylabel("$lnK_d$", fontsize=10, fontweight="bold")
        ax.legend()
        ax.set_title("Vant Hoff") 
        ax.grid(ls=":")
        fig.tight_layout()
        
        
class IsostericHeat(object):
    
    def __init__(self, T1=273.15, T2=293.15):
        """Class for assessing isosteric heat of adsorption using (Clausius–Clapeyron).            
        Parameters
        ----------
        T1 : float or integer, optional
             Adsorption temperature [K] for the 1st adsorption test, by default 273.15
        T2 : float or integer, optional
             Adsorption temperature [K] for the 2nd adsorption test, by default 293.15
        """   
        self.T1 = T1
        self.T2 = T2
        
    
    def set_inlet(self, x1=x_iheat1, y1=y_iheat1, x2=x_iheat2, y2=y_iheat2):
        """Set inlet parameters for Clausius-Clapeyron equation.
        Parameters
        ----------
        x1  : 1d array of floats, optional
              Absolute pressure at kPa for adsorption isotherm at T1
        y1  : 1d array of floats, optional
              Equilibrium adsorption capacity [mmol/g] for adsorption isotherm at T1
        x2  : 1d array of floats, optional
              Absolute pressure at kPa for adsorption isotherm at T2
        y1  : 1d array of floats, optional
              Equilibrium adsorption capacity [mmol/g] for adsorption isotherm at T2    
        """
        x1 = np.array(x1)
        y1 = np.array(y1)
        x2 = np.array(x2)
        y2 = np.array(y2)
        
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        
        if y1.max()>y2.max():
            self.nfit = np.linspace(min(y2)*10, max(y2))
        else:
            self.nfit = np.linspace(min(y1)*10, max(y1))
    
    def freundlich_langmuir(self, x, a, b, c):
        x = np.array(x)
        return a*b*x**c/(1+b*x**c)

    def fl_params_at_T1(self):
        FitParams, FitCov = curve_fit(self.freundlich_langmuir,
                                      self.x1, 
                                      self.y1,
                                      np.array([1, 1, 1]),
                                      bounds=(0, [np.inf, np.inf, np.inf]))
        self.a = FitParams[0]
        self.b = FitParams[1]
        self.c = FitParams[2]
        return {"a at T1 [mmol/g]": self.a,
                "b at T1 [1/kPa^c]": self.b,
                "c at T1 [dimensionless]": self.c,
                }
    
    def fl_curve_at_T1(self, x):
        yfit = self.freundlich_langmuir(x, 
                             self.fl_params_at_T1()["a at T1 [mmol/g]"],
                             self.fl_params_at_T1()["b at T1 [1/kPa^c]"],
                             self.fl_params_at_T1()["c at T1 [dimensionless]"])
        return yfit

    def fl_params_at_T2(self):
        FitParams, FitCov = curve_fit(self.freundlich_langmuir,
                                      self.x2, 
                                      self.y2,
                                      np.array([1, 1, 1]),
                                      bounds=(0, [np.inf, np.inf, np.inf]))
        self.a = FitParams[0]
        self.b = FitParams[1]
        self.c = FitParams[2]
        return {"a at T2 [mmol/g]": self.a,
                "b at T2 [1/kPa^c]": self.b,
                "c at T2 [dimensionless]": self.c,
                }
    
    def fl_curve_at_T2(self, x):
        yfit = self.freundlich_langmuir(x, 
                             self.fl_params_at_T2()["a at T2 [mmol/g]"],
                             self.fl_params_at_T2()["b at T2 [1/kPa^c]"],
                             self.fl_params_at_T2()["c at T2 [dimensionless]"])
        return yfit
    
        
    def plot_freundlich_langmuir(self):
        plt.figure(dpi = 200)
        yfit1 = self.fl_curve_at_T1(self.x1)
        yfit2 = self.fl_curve_at_T2(self.x2)
        plt.plot(self.x1, self.y1, "ko", mfc="none", label = "$T_1={}K$".format(self.T1))
        plt.plot(self.x1, yfit1, "k:", label = "$FL~model$")
        plt.plot(self.x2, self.y2, "ro", mfc="none", label = "$T_2={}K$".format(self.T2))
        plt.plot(self.x2, yfit2, "r:", label = "$FL~model$")
        plt.xlabel("$Absolute~pressure~p~[kPa]$", fontsize=12)
        plt.ylabel("$Amount~adsorbed~n~[mmol/g]$", fontsize=12)
        plt.legend()
        plt.grid(ls=":")  
        plt.tight_layout()
        
    def assess_fit(self):
        yfit1 = self.fl_curve_at_T1(self.x1)
        yfit2 = self.fl_curve_at_T2(self.x2)
        y_observed1 = self.y1
        y_observed2 = self.y2
        n1 = len(y_observed1)
        n2 = len(y_observed2)
         
        R_fl_T1 = 1-((np.sum((y_observed1 - yfit1)**2))/\
                   (np.sum((y_observed1 - np.mean(y_observed1))**2)))\
                   *((n1-1)/(n1-len(self.fl_params_at_T1().items())))
        
        R_fl_T2 = 1-((np.sum((y_observed2 - yfit2)**2))/\
                   (np.sum((y_observed2 - np.mean(y_observed2))**2)))\
                   *((n2-1)/(n2-len(self.fl_params_at_T2().items())))
        
        return {"Freundlich-Langmuir R2 for T1": R_fl_T1, 
                "Freundlich-Langmuir R2 for T2": R_fl_T2}
    
    def isosteric(self, n, a, b, c):
        n = np.array(n)
        return (n/(a*b-n*b))**(1/c)
       
    def clausius_clapeyron(self):
        params_T1 = np.array(list(self.fl_params_at_T1().values()))
        params_T2 = np.array(list(self.fl_params_at_T2().values())) 
        pkPa_T1 = self.isosteric(self.nfit, *params_T1)
        pkPa_T2 = self.isosteric(self.nfit, *params_T2)
        iso_heat = -8.314*np.log(pkPa_T2/pkPa_T1)*self.T1*self.T2/(self.T2-self.T1)
        return -iso_heat/1000
    
    def plot_isoHeat_vs_mmol(self):
        plt.figure(dpi = 200)
        yfit = self.clausius_clapeyron()
        plt.plot(self.nfit, yfit, 'c.', label = "Isosteric heat via Freundlich-Langmuir fit")
        plt.ylabel("$-ΔH_{ads}~[kJ/mol]$", fontsize=12)
        plt.xlabel("$Amount~adsorbed~n~[mmol/g]$", fontsize=12)
        plt.tight_layout()
        plt.legend()
        plt.grid(ls=":")
        plt.ylim(0, yfit.max())
        if self.y1.max()>self.y2.max():
            plt.xlim(0, self.y1.max())
        else:
            plt.xlim(0, self.y2.max())
            
    def all_params(self):
        
        def get_params(*args):
            params = np.vstack(args)
            return params
        
        all_p = get_params(list(self.fl_params_at_T1().items()), 
                           list(self.fl_params_at_T2().items()))
        
        df = pd.DataFrame(all_p, columns = ["Parameters", "Values"])
        return df
            
    def get_dataframe(self):
        params_T1 = np.array(list(self.fl_params_at_T1().values()))
        params_T2 = np.array(list(self.fl_params_at_T2().values()))
        pkPa_T1 = self.isosteric(self.nfit, *params_T1)
        pkPa_T2 = self.isosteric(self.nfit, *params_T2)
        yfit = self.clausius_clapeyron()
        df = pd.DataFrame()
        df["n [mmol/g]"] = self.nfit
        df["pkPa_T1"] = pkPa_T1
        df["pkPa_T2"] = pkPa_T2
        df["ΔH ads [kJ/mol]"] = yfit
        return df
    
    def to_excel(self, filename, **options):
        """Saves the pandas.DataFrame of profiles in an Excel file.
        Parameters
        ----------
        filename : str
            Name of destination file without suffix .xlsx.
        """
        path = filename + '.xlsx'
        with pd.ExcelWriter(path) as writer:
            self.get_dataframe().to_excel(writer, sheet_name="IsoHeat")


class Adsorbent_ScaleUp(object):
  
    def __init__(self, Mr=34.1, y=0.001):
        """Class for calculating the adsorbent quantity in a scaled-up adsorption bed.
        Parameters
        ----------
        Mr : float, optional
             Molar mass of adsorbed molecule [g/mol], by default 34.1 for H2S
        y  : float, optional
             Molar fraction [dimensionless], by default 0.001
        """
        self.Mr = Mr
        self.y = y

    
    def exp_unit(self, ads_capacity=125, exp_ads_time=500):
        """
        Parameters
        ----------
        ads_capacity  : float, optional
                        Adsorption capacity of target molecule (TM) [mg_TM/g_ads], by default 125 mg_TM/g_ads
        exp_ads_time  : float, optional
                        Adsorption time in experimental unit, by default 500 min
        """
        self.ads_capacity = ads_capacity
        self.exp_ads_time = exp_ads_time
        
    def pilot_unit(self, total_flow_rate=2500, pilot_ads_time=2000, print_results = None):
        """
        Parameters
        ----------
        total_flow_rate  : float, optional
                           Total flow rate in pilot unit, by default 2500 mL/min
        pilot_ads_time   : float, optional
                           Desired adsorption time in pilot unit, by default 2000 min
        """
        self.total_flow_rate = total_flow_rate
        self.pilot_ads_time = pilot_ads_time
        self.pilot_unit_target_molecule_flow_rate = self.total_flow_rate*self.y
        self.pilot_unit_target_molecule_flow_rate = self.pilot_unit_target_molecule_flow_rate*(self.Mr/22.4)
        self.mg_adsorbed_quantity_target_molecule = self.pilot_unit_target_molecule_flow_rate*self.pilot_ads_time
        self.quantity_ads_in_pilot_unit = self.mg_adsorbed_quantity_target_molecule/self.ads_capacity

        print("The quantity of adsorbent in the pilot unit in (g) is:", self.quantity_ads_in_pilot_unit)
