import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

class Isotherms(object):
    
    def __init__(self, x, y, P = 1, Mr = 44.01, T = 298.15, R = 8.205e-5):
        """Class for 6 different adsorption isotherms.
        Parameters
        ----------
        x  : 1d array of floats
             Equilibrium concentration dimensionless [%]
        y  : 1d array of floats
             Equilibrium adsorption capacity [mg/g]
        P  : float or integer, optional
             Adsorption pressure [atm], by default 1
        Mr : float, optional
             Molar mass of adsorbed molecule [g/mol], by default 44.01 for CO2
        T  : float, optional
             Adsorption temperature [K], by default 298.15
        R  : float, optional
             Gas constant [atm.m3/mol/K], by default 8.205e-5
        
        Calculated Parameters
        ----------
        factor  : 1d array of floats
                  Factor to convert [%] to [mg/L]
        x_obs   : 1d array of floats
                  Converted equilibrium concentration [mg/L]
        xfit    : 1d array of floats
                  Predictor for the models tested [mg/L]    
        """
        self.x = x
        self.y = y
        self.P = P
        self.Mr = Mr
        self.T = T
        self.R = R
        self.factor = (self.P*self.Mr)/self.R/self.T
        self.x_obs = self.factor*self.x
        self.xfit = np.linspace(min(self.x_obs), max(self.x_obs), 50)
    
    def langmuir(self, x, k, q):
        x = np.array(x)
        return k*q*x/(1+k*x)

    def langmuir_params(self):
        FitParams, FitCov = curve_fit(self.langmuir,
                                      self.x_obs, 
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.q = FitParams[1]
        return {'K_Langmuir [L/mg]': self.k,
                'qmax_Langmuir [mg/g]': self.q}
    
    def langmuir_curve(self, x):
        yfit = self.langmuir(x, 
                             self.langmuir_params()['K_Langmuir [L/mg]'],
                             self.langmuir_params()['qmax_Langmuir [mg/g]'])
        return yfit
    
    def plot_langmuir_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.langmuir_curve(self.xfit), 'k--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('Langmuir fit')
        ax.grid(ls=":")
    
    
    def freundlich(self, x, k, n):
        x = np.array(x)
        return k*x**(1/n)
    
    def freundlich_params(self):
        FitParams, FitCov = curve_fit(self.freundlich,
                                      self.x_obs, 
                                      self.y,
                                      np.array([1, 1]), 
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.n = FitParams[1]
        return {'K_Freundlich [L^(1/n_Freundlich)·mg^(1-1/n_Freundlich)·g^-1]': 
                self.k, 'n_Freundlich': self.n}
    
    def freundlich_curve(self, x):
        yfit = self.freundlich(x, 
                               self.freundlich_params()['K_Freundlich [L^(1/n_Freundlich)·mg^(1-1/n_Freundlich)·g^-1]'],
                               self.freundlich_params()['n_Freundlich'])
        return yfit
    
    def plot_freundlich_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.freundlich_curve(self.xfit), 'k--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('Freundlich fit')
        ax.grid(ls=":")
    
    
    def temkin(self, x, a, b):
        x = np.array(x)
        return (8.314*298.15/b)*(np.log(a*x))
    
    def temkin_params(self):
        FitParams, FitCov = curve_fit(self.temkin,
                                      self.x_obs, 
                                      self.y,
                                      np.array([1, 1]), 
                                      bounds=(0, [np.inf, np.inf]))
        self.a = FitParams[0]
        self.b = FitParams[1]
        return {'A_Temkin [L/mg]': self.a,
                'B_Temkin [J/mol]': self.b}
    
    def temkin_curve(self, x):
        yfit = self.temkin(x, 
                           self.temkin_params()['A_Temkin [L/mg]'],
                           self.temkin_params()['B_Temkin [J/mol]'])
        return yfit
    
    def plot_temkin_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.temkin_curve(self.xfit), 'k--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('Temkin fit')
        ax.grid(ls=":")
    
    
    def toth(self, x, k, n, q):
        x = np.array(x)
        return q*x/(k+x**n)**(1/n)
    
    def toth_params(self):
        FitParams, FitCov = curve_fit(self.toth,
                                      self.x_obs, 
                                      self.y,
                                      np.array([1, 1, 1]), 
                                      bounds=(0, [np.inf, np.inf, np.inf]))
        self.k = FitParams[0]
        self.n = FitParams[1]
        self.q = FitParams[2]
        return {'K_Toth [mg^z·L^-n_toth]': self.k,
                'n_Toth': self.n,
                'qmax_Toth [mg/g]': self.q}
    
    def toth_curve(self, x):
        yfit = self.toth(x, 
                         self.toth_params()['K_Toth [mg^z·L^-n_toth]'],
                         self.toth_params()['n_Toth'],
                         self.toth_params()['qmax_Toth [mg/g]'])
        return yfit

    def plot_toth_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.toth_curve(self.xfit), 'k--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('Toth fit')
        ax.grid(ls=":")
        
    
    def sips(self, x, k, n, q):
        x = np.array(x)
        return (q*k*x**(1/n))/(1+k*x**(1/n))
    
    def sips_params(self):
        FitParams, FitCov = curve_fit(self.sips, 
                                      self.x_obs, 
                                      self.y,
                                      np.array([1, 1, 1]), 
                                      bounds=(0, [np.inf, np.inf, np.inf]))
        self.k = FitParams[0]
        self.n = FitParams[1]
        self.q = FitParams[2]
        return {'K_Sips [L^n_Sips·mg^-n_Sips]': self.k,
                'n_Sips': self.n,
                'qmax_Sips [mg/g]': self.q}

    def sips_curve(self, x):
        yfit = self.sips(x, 
                         self.sips_params()['K_Sips [L^n_Sips·mg^-n_Sips]'],
                         self.sips_params()['n_Sips'],
                         self.sips_params()['qmax_Sips [mg/g]'])
        return yfit   
    
    def plot_sips_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.sips_curve(self.xfit), 'k--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('Sips fit')  
        ax.grid(ls=":")
        
    
    def dubinin_radushkevich(self, x, E, q):
        x = np.array(x)
        return q*np.exp((8.314*298.15*np.log(1+1/x))**2/(-2*E**2))

    def dubinin_radushkevich_params(self):
        FitParams, FitCov = curve_fit(self.dubinin_radushkevich, 
                                      self.x_obs, 
                                      self.y,
                                      np.array([1, 1]), 
                                      bounds=(0, [np.inf, np.inf]))
        self.E = FitParams[0]
        self.q = FitParams[1]
        return {'E_DR [J/mol]': self.E, 
                'qmax_DR [mg/g]': self.q}
    
    def dr_curve(self, x):
        yfit = self.dubinin_radushkevich(x, 
                                         self.dubinin_radushkevich_params()['E_DR [J/mol]'],
                                         self.dubinin_radushkevich_params()['qmax_DR [mg/g]'])
        return yfit

    def plot_dr_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.dr_curve(self.xfit), 'k--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('DR fit') 
        ax.grid(ls=":")
    

    def plot_all_models(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.langmuir_curve(self.xfit), 'r--', mfc = 'none', label = 'Langmuir')
        ax.plot(self.xfit, self.freundlich_curve(self.xfit), 'b--', mfc = 'none', label = 'Freundlich')
        ax.plot(self.xfit, self.temkin_curve(self.xfit), 'g--', mfc = 'none', label = 'Temkin')
        ax.plot(self.xfit, self.toth_curve(self.xfit), 'y--', mfc = 'none', label = 'Toth')
        ax.plot(self.xfit, self.sips_curve(self.xfit), 'm--', mfc = 'none', label = 'Sips')
        ax.plot(self.xfit, self.dr_curve(self.xfit), 'c--', mfc = 'none', label = 'DR')
        ax.set_xlabel("C$_e$ [mg/L]", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("q$_e$ [mg/g]", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('All models') 
        ax.grid(ls=":")        


    def assess_fit(self):
        y_observed = self.y
        y_langmuir = self.langmuir_curve(self.x_obs)
        y_freundlich = self.freundlich_curve(self.x_obs)
        y_temkin = self.temkin_curve(self.x_obs)
        y_toth = self.toth_curve(self.x_obs)
        y_sips = self.sips_curve(self.x_obs)
        y_DR = self.dr_curve(self.x_obs)
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_langmuir)
        R_langmuir = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_freundlich)
        R_freundlich = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_temkin)
        R_temkin = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_toth)
        R_toth = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_sips)
        R_sips = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_DR)
        R_dr = r**2   
        return {"Langmuir R2": R_langmuir, 
                "Freundlich R2": R_freundlich,
                "Temkin R2": R_temkin,
                "Toth R2": R_toth,
                "Sips R2": R_sips,
                "DR R2": R_dr}
    
    
    def best_fit(self):
        model = max(self.assess_fit(), key=self.assess_fit().get)
        value = self.assess_fit().get(model)
        return print("The best model is that of", model, "=", value)
    
    
    def all_params(self):
        
        def get_params(*args):
            params = np.vstack(args)
            return params
        
        all_p = get_params(list(self.langmuir_params().items()), 
                           list(self.freundlich_params().items()),
                           list(self.temkin_params().items()),
                           list(self.toth_params().items()),
                           list(self.sips_params().items()),
                           list(self.dubinin_radushkevich_params().items()))
        
        df = pd.DataFrame(all_p, columns = ['Parameters', 'Values'])
        return df
    


class Kinetics(object):

    def __init__(self, x, y):
        """Class for 6 different adsorption kinetic models.
        Parameters
        ----------
        x  : 1d array of floats
             Adsorption time [min]
        y  : 1d array of floats
             Accumulative adsorption capacity [mg/g]   
        """
        self.x = x
        self.y = y
        
    def pfo(self, x, k, q):
        x = np.array(x)
        return q*(1-np.exp(-k*x))
    
    def pfo_params(self):
        FitParams, FitCov = curve_fit(self.pfo,
                                      self.x,
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.q = FitParams[1]
        return {'k_pfo [1/min]': self.k,
                'qmax_PFO [mg/g]': self.q}
    
    def pfo_curve(self, x):
        yfit = self.pfo(x, 
                        self.pfo_params()['k_pfo [1/min]'],
                        self.pfo_params()['qmax_PFO [mg/g]'])
        return yfit
    
    def plot_pfo_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.pfo_curve(self.x), 'r--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('PFO fit')
        ax.grid(ls=":")
     

    def pso(self, x, k, q):
        x = np.array(x)
        return k*q**2*x/(1+k*q*x)

    def pso_params(self):
        FitParams, FitCov = curve_fit(self.pso,
                                      self.x,
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.q = FitParams[1]
        return {'k_pso [g mg^-1 min^-1]': self.k,
                'qmax_PSO [mg/g]': self.q}
    
    def pso_curve(self, x):
        yfit = self.pso(x, 
                        self.pso_params()['k_pso [g mg^-1 min^-1]'],
                        self.pso_params()['qmax_PSO [mg/g]'])
        return yfit

    def plot_pso_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.pso_curve(self.x), 'r--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('PSO fit')
        ax.grid(ls=":")

        
    def weber_morris(self, x, k, c):
        x = np.array(x)
        return k*x**(0.5)+c
    
    def weber_morris_params(self):
        FitParams, FitCov = curve_fit(self.weber_morris,
                                      self.x, 
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.c = FitParams[1]
        return {'k_wm [mg g^-1 min^-0.5]': self.k,
                'C': self.c}
    
    def weber_morris_curve(self, x):
        yfit = self.weber_morris(x, 
                                 self.weber_morris_params()['k_wm [mg g^-1 min^-0.5]'],
                                 self.weber_morris_params()['C'])
        return yfit

    def plot_weber_morris_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.weber_morris_curve(self.x), 'r--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Weber-Morris fit')
        ax.grid(ls=":")
     

    def avrami(self, x, k, q, n):
        x = np.array(x)
        return q*(1-np.exp(-k*x))**n

    def avrami_params(self):
        FitParams, FitCov = curve_fit(self.avrami,
                                      self.x, 
                                      self.y,
                                      np.array([1, 1, 1]),
                                      bounds=(0, [np.inf, np.inf, np.inf]))
        self.k = FitParams[0]
        self.q = FitParams[1]
        self.n = FitParams[2]
        return {'k_avrami [1/min]': self.k,
                'qmax_avrami [mg/g]': self.q,
                'n_avrami': self.n}

    def avrami_curve(self, x):
        yfit = self.avrami(x, 
                           self.avrami_params()['k_avrami [1/min]'],
                           self.avrami_params()['qmax_avrami [mg/g]'],
                           self.avrami_params()['n_avrami'])
        return yfit

    def plot_avrami_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.avrami_curve(self.x), 'r--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Avrami fit')
        ax.grid(ls=":")
        
    
    def bangham(self, x, k, q, n):
        x = np.array(x)
        return q*(1-np.exp(-k*x**n))
    
    def bangham_params(self):
        FitParams, FitCov = curve_fit(self.bangham,
                                      self.x, 
                                      self.y,
                                      np.array([1, 1, 1]),
                                      bounds=(0, [np.inf, np.inf, np.inf]))
        self.k = FitParams[0]
        self.q = FitParams[1]
        self.n = FitParams[2]
        return {'k_bangham [1/min^n]': self.k,
                'qmax__bangham [mg/g]': self.q,
                'n_bangham': self.n}
    
    def bangham_curve(self, x):
        yfit = self.bangham(x, 
                           self.bangham_params()['k_bangham [1/min^n]'],
                           self.bangham_params()['qmax__bangham [mg/g]'],
                           self.bangham_params()['n_bangham'])
        return yfit

    def plot_bangham_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.bangham_curve(self.x), 'r--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Bangham fit')
        ax.grid(ls=":")

        
    def elovich(self, x, a, b):
        x = np.array(x)
        return (1/b)*np.log(a*b*x)

    def elovich_params(self):
        FitParams, FitCov = curve_fit(self.elovich,
                                      self.x[1:], 
                                      self.y[1:],
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.a = FitParams[0]
        self.b = FitParams[1]
        return {'a_elovich [mg/g/min]': self.a,
                'b__elovich [g/mg]': self.b}
    
    def elovich_curve(self, x):
        yfit = self.elovich(x, 
                           self.elovich_params()['a_elovich [mg/g/min]'],
                           self.elovich_params()['b__elovich [g/mg]'])
        return yfit

    def plot_elovich_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.elovich_curve(self.x), 'r--', mfc = 'none', label = 'Predicted')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Elovich fit')
        ax.set_ylim(0, np.max(self.y)+20)
        ax.grid(ls=":")
        
    def plot_all_models(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.pfo_curve(self.x), 'r--', mfc = 'none', label = 'PFO')
        ax.plot(self.x, self.pso_curve(self.x), 'b--', mfc = 'none', label = 'PSO')     
        ax.plot(self.x, self.weber_morris_curve(self.x), 'g--', mfc = 'none', label = 'Weber-Morris')
        ax.plot(self.x, self.avrami_curve(self.x), 'y--', mfc = 'none', label = 'Avrami')
        ax.plot(self.x, self.bangham_curve(self.x), 'm--', mfc = 'none', label = 'Bangham')
        ax.plot(self.x, self.elovich_curve(self.x), 'c--', mfc = 'none', label = 'Elovich')
        ax.set_xlabel("Time [min]", fontsize=10, fontweight='bold')
        ax.set_ylabel("q$_t$ [mg/g]", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('All models') 
        ax.set_ylim(0, np.max(self.y)+20)
        ax.grid(ls=":")


    def assess_fit(self):
        y_observed = self.y
        y_pfo = self.pfo_curve(self.x)
        y_pso = self.pso_curve(self.x)
        y_wm = self.weber_morris_curve(self.x)
        y_avrami = self.avrami_curve(self.x)
        y_bangham = self.bangham_curve(self.x)
        y_elovich = self.elovich_curve(self.x)
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_pfo)
        R_pfo = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_pso)
        R_pso = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_wm)
        R_wm = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_avrami)
        R_avrami = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_bangham)
        R_bangham = r**2
        slope, intercept, r, p, std_err = stats.linregress(y_observed, y_elovich)
        R_elovich = r**2  
        return {"PFO R2": R_pfo, 
                "PSO R2": R_pso,
                "WEBER-MORRIS R2": R_wm,
                "AVRAMI R2": R_avrami,
                "BANGHAM R2": R_bangham,
                "ELOVICH R2": R_elovich}

    
    def best_fit(self):
        model = max(self.assess_fit(), key=self.assess_fit().get)
        value = self.assess_fit().get(model)
        return print("The best model is that of", model, "=", value)
    
    
    def all_params(self):
        
        def get_params(*args):
            params = np.vstack(args)
            return params
        
        all_p = get_params(list(self.pfo_params().items()), 
                           list(self.pso_params().items()),
                           list(self.weber_morris_params().items()),
                           list(self.avrami_params().items()),
                           list(self.bangham_params().items()),
                           list(self.elovich_params().items()))
        
        df = pd.DataFrame(all_p, columns = ['Parameters', 'Values'])
        return df

    
class ModifiedArrhenius(object):
    
    def __init__(self, x, y):
        """Class for calculatin Arrhenius parameters.
        Parameters
        ----------
        x  : 1d array of floats
             Adsorption temperatures [K]
        y  : 1d array of floats
             Kinetic constants of the best fitted model at different temperatures   
        """
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
        return {'Ea [J/mol]': self.Ea,
                'A [same as k]': self.A}

    def arrhenius_curve(self, x):
        yfit = self.arrhenius(x, 
                              self.arrhenius_params()['Ea [J/mol]'],
                              self.arrhenius_params()['A [same as k]'])
        return yfit
        
    def plot_arrhenius_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        xfit = np.linspace(min(self.x), max(self.x), 200)
        ax.plot(self.x, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(xfit, self.arrhenius_curve(xfit), 'r--', mfc = 'none', label = 'Arrhenius')
        ax.set_xlabel("Temperature [K]", fontsize=10, fontweight='bold')
        ax.set_ylabel("k", fontsize=10, fontweight='bold')
        ax.set_title("Modified Arrhenius", fontsize=10, fontweight='bold')
        ax.grid(linestyle=':')