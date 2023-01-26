import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from adsorption.ads_data import (x_iso, y_iso, x_kin, y_kin, x_arrh, y_arrh, 
                      x_dyn, y_dyn, x_h, y_h)

class Isotherms(object):
    
    def __init__(self, P = 1, Mr = 44.01, T = 298.15, R = 8.205e-5):
        """Class for evaluating 6 different adsorption isotherms.
        Parameters
        ----------
             Equilibrium adsorption capacity [mg/g]
        P  : float or integer, optional
             Adsorption pressure [atm], by default 1
        Mr : float, optional
             Molar mass of adsorbed molecule [g/mol], by default 44.01 for CO2
        T  : float, optional
             Adsorption temperature [K], by default 298.15
        R  : float, optional
             Gas constant [atm.m3/mol/K], by default 8.205e-5
        """
        self.P = P
        self.Mr = Mr
        self.T = T
        self.R = R
        
    def set_inlet(self, x=x_iso, y=y_iso):
        """Set inlet parameters for isotherm equations.
        Parameters
        ----------
        x  : 1d array of floats, optional
             Equilibrium concentration dimensionless [%]
        y  : 1d array of floats, optional
             Equilibrium adsorption capacity [mg/g]

        Calculated Parameters
        ----------
        factor  : 1d array of floats
                  Factor to convert [%] to [mg/L]
        x_obs   : 1d array of floats
                  Converted equilibrium concentration [mg/L]
        xfit    : 1d array of floats
                  Predictor for the models tested [mg/L]    
        """
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y
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
        ax.plot(self.xfit, self.langmuir_curve(self.xfit), 'k--', label = 'Predicted')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
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
        ax.plot(self.xfit, self.freundlich_curve(self.xfit), 'k--', label = 'Predicted')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
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
        ax.plot(self.xfit, self.temkin_curve(self.xfit), 'k--', label = 'Predicted')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
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
        ax.plot(self.xfit, self.toth_curve(self.xfit), 'k--', label = 'Predicted')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
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
        ax.plot(self.xfit, self.sips_curve(self.xfit), 'k--', label = 'Predicted')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
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
        ax.plot(self.xfit, self.dr_curve(self.xfit), 'k--', label = 'Predicted')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
        ax.legend()
        ax.set_title('DR fit') 
        ax.grid(ls=":")
    
    def plot_all_models(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_obs, self.y, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.xfit, self.langmuir_curve(self.xfit), 'r--', label = 'Langmuir')
        ax.plot(self.xfit, self.freundlich_curve(self.xfit), 'b--', label = 'Freundlich')
        ax.plot(self.xfit, self.temkin_curve(self.xfit), 'g--', label = 'Temkin')
        ax.plot(self.xfit, self.toth_curve(self.xfit), 'y--', label = 'Toth')
        ax.plot(self.xfit, self.sips_curve(self.xfit), 'm--', label = 'Sips')
        ax.plot(self.xfit, self.dr_curve(self.xfit), 'c--', label = 'DR')
        ax.set_xlabel("$C_e$ $[mg/L]$", fontsize = 10, fontweight = 'bold')
        ax.set_ylabel("$q_e$ $[mg/g]$", fontsize = 10, fontweight = 'bold')
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
        n = len(y_observed)
        
        R_langmuir = 1-((np.sum((y_observed - y_langmuir)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.langmuir_params().items())))
        
        R_freundlich = 1-((np.sum((y_observed - y_freundlich)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.freundlich_params().items())))
        
        R_temkin = 1-((np.sum((y_observed - y_temkin)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.temkin_params().items())))
        
        R_toth = 1-((np.sum((y_observed - y_toth)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.toth_params().items())))
        
        
        R_sips = 1-((np.sum((y_observed - y_sips)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.sips_params().items())))
                   
        R_dr = 1-((np.sum((y_observed - y_DR)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.dubinin_radushkevich_params().items())))
        
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
    
    def to_excel(self, filename, **options):
        """Saves the pandas.DataFrame of profiles in an Excel file.
        Parameters
        ----------
        filename : str
            Name of destination file without suffix .xlsx.
        """
        path = filename + '.xlsx'
        with pd.ExcelWriter(path) as writer:
            self.all_params().to_excel(writer, sheet_name='Isotherms')
    

class Kinetics(object):

    def __init__(self):
        """Class for evaluating 6 different adsorption kinetic models.  
        """
        
    def set_inlet(self, x=x_kin, y=y_kin):
        """Set inlet parameters for kinetic equations.
        Parameters
        ----------
        x  : 1d array of floats, optional
             Adsorption time [min]
        y  : 1d array of floats, optional
             Accumulative adsorption capacity [mg/g]   
        """
        x = np.array(x)
        y = np.array(y)
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
        ax.plot(self.x, self.pfo_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
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
        ax.plot(self.x, self.pso_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
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
        ax.plot(self.x, self.weber_morris_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
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
        ax.plot(self.x, self.avrami_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
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
                'qmax_bangham [mg/g]': self.q,
                'n_bangham': self.n}
    
    def bangham_curve(self, x):
        yfit = self.bangham(x, 
                           self.bangham_params()['k_bangham [1/min^n]'],
                           self.bangham_params()['qmax_bangham [mg/g]'],
                           self.bangham_params()['n_bangham'])
        return yfit

    def plot_bangham_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.bangham_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
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
                'b_elovich [g/mg]': self.b}
    
    def elovich_curve(self, x):
        yfit = self.elovich(x, 
                           self.elovich_params()['a_elovich [mg/g/min]'],
                           self.elovich_params()['b_elovich [g/mg]'])
        return yfit

    def plot_elovich_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.elovich_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Elovich fit')
        ax.set_ylim(0, np.max(self.y)+20)
        ax.grid(ls=":")
        
    def plot_all_models(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.pfo_curve(self.x), 'r--', label = 'PFO')
        ax.plot(self.x, self.pso_curve(self.x), 'b--', label = 'PSO')     
        ax.plot(self.x, self.weber_morris_curve(self.x), 'g--', label = 'Weber-Morris')
        ax.plot(self.x, self.avrami_curve(self.x), 'y--', label = 'Avrami')
        ax.plot(self.x, self.bangham_curve(self.x), 'm--', label = 'Bangham')
        ax.plot(self.x, self.elovich_curve(self.x), 'c--', label = 'Elovich')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$q_t$ $[mg/g]$", fontsize=10, fontweight='bold')
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
        y_elovich = self.elovich_curve(self.x[1:])
        n = len(y_observed)
               
        R_pfo = 1-((np.sum((y_observed - y_pfo)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.pfo_params().items())))
        
        R_pso = 1-((np.sum((y_observed - y_pso)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.pso_params().items())))
        
        R_wm = 1-((np.sum((y_observed - y_wm)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.weber_morris_params().items())))
                   
        R_avrami = 1-((np.sum((y_observed - y_avrami)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.avrami_params().items())))
        
        R_bangham = 1-((np.sum((y_observed - y_bangham)**2))/\
                   (np.sum((y_observed - np.mean(y_observed))**2)))\
                   *((n-1)/(n-len(self.bangham_params().items())))
        
        
        R_elovich = 1-((np.sum((y_observed[1:] - y_elovich)**2))/\
                   (np.sum((y_observed[1:] - np.mean(y_observed[1:]))**2)))\
                   *((n-1)/(n-len(self.elovich_params().items())))  
        
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

    def to_excel(self, filename, **options):
        """Saves the pandas.DataFrame of profiles in an Excel file.
        Parameters
        ----------
        filename : str
            Name of destination file without suffix .xlsx.
        """
        path = filename + '.xlsx'
        with pd.ExcelWriter(path) as writer:
            self.all_params().to_excel(writer, sheet_name='Kinetics')

    
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
        ax.plot(xfit, self.arrhenius_curve(xfit), 'g--', label = 'Arrhenius')
        ax.set_xlabel("$Temperature$ $[K]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$k$", fontsize=10, fontweight='bold')
        ax.set_title("Modified Arrhenius", fontsize=10)
        ax.grid(linestyle=':')
        
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
        """Class for calculating 3 different empirical adsorpion dynamic models.            
        Parameters
        ----------
        C  : float or integer, optional
             Initial concentration of adsorbed molecule [%], by default 0.1 (10% CO2)
        Mr : float, optional
             Molar mass of adsorbed molecule [g/mol], by default 44.01 for CO2
        T  : float or integer, optional
             Adsorption temperature [K], by default 298.15
        P  : float or integer, optional
             Adsorption pressure [atm], by default 298.15
        h  : float or integer, optional
             Height of adsorption bed [cm], by default 2
        r  : float or integer, optional
             Adsorption bed inner radius [cm], by default 0.9
        Q  : integer, optional
             Total flow rate of gas mixture [ml/min], by default 0.9
        W  : float or integer, optional
             Weight of adsorbent used [g], by default 1
        U  : float or integer, optional
             Considering part of adsorption process, by default 0.1 (10% of Ct/C0 for Adams-Bohart eq.)
        R  : float, optional
             Gas constant [atm.m3/mol/K], by default 8.205e-5
        
        Calculated Parameters
        ----------
        c0      : float
                  Converted equilibrium concentration [mg/ml]
        A       : float
                  Cross-sectional area of the reactor bed [cm2]
        v       : float
                  Superficial velocity [cm/min]           
        """             
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
        self.c0 = self.C*(self.P*self.Mr)/self.R/self.T/1000
        self.A = np.pi*(self.r**2)
        self.v = self.Q/self.A

    def set_inlet(self, x=x_dyn, y=y_dyn):
        """Set inlet parameters for dynamic equations.
        Parameters
        ----------
        x  : 1d array of floats, optional
             Adsorption time [min]
        y  : 1d array of floats, optional
             Ct/C0 dimensionless concentration 
        
        Calculated Parameters
        ----------
        yy      : 1d array of floats
                  Initial stage of adsorption process for Adams-Bohart equation (Ct/C0 = 0.1)
        xx      : 1d array of floats
                  Minutes correspond to Ct/C0 = 0.1   
        """
        x = np.array(x)
        y = np.array(y)
        self.x = x
        self.y = y
        self.yy = self.y[self.y <= self.U]
        self.xx = self.x[np.where(self.y <= self.U)[0]]

    def thomas(self, x, k, q):
        x = np.array(x)
        return 1/(1 + np.exp(k*q*(self.W/self.Q) - k*self.c0*x)) 

    def thomas_params(self):
        FitParams, FitCov = curve_fit(self.thomas,
                                      self.x, 
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.q = FitParams[1]
        return {'k_thomas [ml/mg/min]': self.k,
                'qmax_thomas [mg/g]': self.q}       

    def thomas_curve(self, x):
        yfit = self.thomas(x, 
                           self.thomas_params()['k_thomas [ml/mg/min]'],
                           self.thomas_params()['qmax_thomas [mg/g]'])
        return yfit

    def plot_thomas_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.thomas_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$C_t/C_0$", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Thomas fit')
        ax.grid(ls=":")
        
    def yoon_nelson(self, x, k, tau):
        x = np.array(x)
        return np.exp(k*(x-tau))/(1+np.exp(k*(x-tau)))         
        
    def yoon_nelson_params(self):
        FitParams, FitCov = curve_fit(self.yoon_nelson,
                                      self.x, 
                                      self.y,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.tau = FitParams[1]
        return {'k_yoon_nelson [1/min]': self.k,
                'tau_yoon_nelson [min]': self.tau} 

    def yoon_nelson_curve(self, x):
        yfit = self.yoon_nelson(x, 
                                self.yoon_nelson_params()['k_yoon_nelson [1/min]'],
                                self.yoon_nelson_params()['tau_yoon_nelson [min]'])
        return yfit

    def plot_yoon_nelson_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.yoon_nelson_curve(self.x), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$C_t/C_0$", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Yoon-Nelson fit')
        ax.grid(ls=":")

    def adams_bohart(self, x, k, N0):
        x = np.array(x)
        return np.exp(k*self.c0*x - k*N0*(self.h/self.v))
    
    def adams_bohart_params(self):
        FitParams, FitCov = curve_fit(self.adams_bohart,
                                      self.xx, 
                                      self.yy,
                                      np.array([1, 1]),
                                      bounds=(0, [np.inf, np.inf]))
        self.k = FitParams[0]
        self.N0 = FitParams[1]
        return {'k_adams_bohart [ml/mg/min]': self.k,
                'N0_adams_bohart [mg/ml]': self.N0}
    
    def adams_bohart_curve(self, x):
        yfit = self.adams_bohart(x, 
                                self.adams_bohart_params()['k_adams_bohart [ml/mg/min]'],
                                self.adams_bohart_params()['N0_adams_bohart [mg/ml]'])
        return yfit

    def plot_adams_bohart_fit(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.xx, self.adams_bohart_curve(self.xx), 'r--', label = 'Predicted')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$C_t/C_0$", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Adams-Bohart fit')
        ax.grid(ls=":")

    def plot_all_models(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x, self.y, 'k:', mfc = 'none', label = 'Observed')
        ax.plot(self.x, self.thomas_curve(self.x), 'r--', label = 'Thomas')
        ax.plot(self.x, self.yoon_nelson_curve(self.x), 'b--', label = 'Yoon-Nelson')
        ax.plot(self.xx, self.adams_bohart_curve(self.xx), 'm--', label = 'Adams-Bohart')
        ax.set_xlabel("$Time$ $[min]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$C_t/C_0$", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('All models') 
        ax.grid(ls=":")
        
    def assess_fit(self):
        y_obs1 = self.y
        y_obs2 = self.yy
        y_thomas = self.thomas_curve(self.x)
        y_yoon_nelson = self.yoon_nelson_curve(self.x)
        y_adams_bohart = self.adams_bohart_curve(self.xx)
        n1 = len(y_obs1)
        n2 = len(y_obs2)
        
        R_thomas = 1-((np.sum((y_obs1 - y_thomas)**2))/\
                   (np.sum((y_obs1 - np.mean(y_obs1))**2)))\
                   *((n1-1)/(n1-len(self.thomas_params().items())))
                   
        R_yoon_nelson = 1-((np.sum((y_obs1 - y_yoon_nelson)**2))/\
                   (np.sum((y_obs1 - np.mean(y_obs1))**2)))\
                   *((n1-1)/(n1-len(self.yoon_nelson_params().items())))
        
        R_adams_bohart = 1-((np.sum((y_obs2 - y_adams_bohart)**2))/\
                   (np.sum((y_obs2 - np.mean(y_obs2))**2)))\
                   *((n2-1)/(n2-len(self.adams_bohart_params().items())))
        
        return {"THOMAS R2": R_thomas, 
                "YOON-NELSON R2": R_yoon_nelson,
                "ADAMS-BOHART R2": R_adams_bohart}
    
    def best_fit(self):
        model = max(self.assess_fit(), key=self.assess_fit().get)
        value = self.assess_fit().get(model)
        return print("The best model is that of", model, "=", value)
    
    def all_params(self):
        
        def get_params(*args):
            params = np.vstack(args)
            return params
        
        all_p = get_params(list(self.thomas_params().items()), 
                           list(self.yoon_nelson_params().items()),
                           list(self.adams_bohart_params().items()))
        
        df = pd.DataFrame(all_p, columns = ['Parameters', 'Values'])
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
            self.all_params().to_excel(writer, sheet_name='AdsDynamics')


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
        """Set inlet parameters for modified Vant Hoff equation.
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
        return {'enthalpy [J/mol]': enthlapy,
                'entropy [J/mol/K]': entropy,
                'R2': r**2, 
                'slope': slope, 
                'intercept': intercept}
        
    def vant_hoff_line(self, x):
        yfit = list(map(lambda x: self.vant_hoff_params()['slope']*x 
                        + self.vant_hoff_params()['intercept'], 1/x))
        return np.array(yfit)
    
    def plot_vant_hoff(self):
        fig, ax = plt.subplots(figsize = (6,4), dpi = 200)
        ax.plot(self.x_inv, self.lnKd, 'ko', mfc = 'none', label = 'Observed')
        ax.plot(self.x_inv, self.vant_hoff_line(self.x), 'r--', label = 'Vant Hoff')
        ax.set_xlabel("$1/T$ $[1/K]$", fontsize=10, fontweight='bold')
        ax.set_ylabel("$lnK_d$", fontsize=10, fontweight='bold')
        ax.legend()
        ax.set_title('Vant Hoff') 
        ax.grid(ls=":")
