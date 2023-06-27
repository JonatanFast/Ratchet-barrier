# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 21:23:08 2022

@author: jo0651fa
"""
import numpy as np
import matplotlib.pyplot  as plt
import warnings
from tqdm import tqdm
import json
import tkinter as tk
from tkinter import filedialog as fd


# 

#############
# Functions #
#############
def get_file(title='Select File', extension=''):
    """
    Opens a user file selection window and returns selected file directory.
    
    Parameters
    ----------
    title : str
        Window title.

    Returns
    -------
    filename : str
        User selected file directory.

    """
    root = tk.Tk() # Define tkinter window
    root.attributes("-topmost", True) #Force file selection dialog on top
    # NOTE: root.attributes is Windows specific!
    root.withdraw() #Hide tkinter window
    # Open file selection diolog; if certaine xtension specified: Check for correct file extension
    if extension == '.mat':
        filename = fd.askopenfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.mat',
                                            filetypes=(("matlab workspace file", "*.mat"),))
    if extension == '.json':
        filename = fd.askopenfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.json',
                                            filetypes=(("json file", "*.json"),))
    if extension == '.txt':
        filename = fd.askopenfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.txt',
                                            filetypes=(("text file", "*.txt"),))
    if not extension:
        filename = fd.askopenfilename(parent=root, title='{}'.format(title))        
    return filename

def get_savefile(title='Save file as', extension=''):
    
    """
    Opens a user file selection window and returns selected save-file directory.
    
    Parameters
    ----------
    title : str
        Window title. Optional.
    extension : str
        Forces file extension in str. So far only '.json' possible.
        
    Returns
    -------
    filename : str
        User selected save-file directory.

    """
    root = tk.Tk() # Define tkinter window
    root.attributes("-topmost", True) #Force file selection dialog on top
    # NOTE: root.attributes is Windows specific!
    root.withdraw() #Hide tkinter window
    # Open file selectiopn diolog, depending on specified file extension:
    if not extension:
        savefilename = fd.asksaveasfilename(parent=root, title='{}'.format(title))
    
    elif extension == '.json':
        savefilename = fd.asksaveasfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.json',
                                            filetypes=(("json file", "*.json"),))
    elif extension == '.svg':
        savefilename = fd.asksaveasfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.svg',
                                            filetypes=(("svg file", "*.svg"),))
    elif extension == '.pdf':
        savefilename = fd.asksaveasfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.pdf',
                                            filetypes=(("pdf file", "*.pdf"),))
    elif extension == 'image':
        savefilename = fd.asksaveasfilename(parent=root, title='{}'.format(title),
                                            defaultextension='.png',
                                            filetypes=(("png file", "*.png"),("svg file", "*.svg")))
    # # If none of specified extensions: Return to default - all extensions.
    # else:
    #     savefilename = fd.asksaveasfilename(parent=root, title='{}'.format(title))
        
    return savefilename

def find_nearest(array, value): #
    '''
    Returns index of the element in array which is closest to specified value 
    Parameters
    ----------
    array : np.array
        
    value : float
        
    Returns
    -------
    idx : int       
    '''    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx #array[idx]

def WKB_model(mu, kDeltaT_G, kDeltaT_S, V_model, Ampl, L, dE = 0.1, T = 77, U_top = 340, m_eff=0.023):
    '''
    Uses the Landauer-Buttiker scattering theory to determine current based on fermidistributionon each side and transmission fucniotn T
    transmission T is determined by WKB model for triangular barrier, and only accounts for tunneling, T=1 at energies above the barrier
    default values:
        L: length of barrier in nm
        dE: Energy step size, in [meV], for Landau-Buttiker integration
        kT: base temperature in meV, temp in K multiplied by Boltzman in meV/K
        U_top: barrier height in meV
        Ampl: scalling factor between experimental and calculated current, determined from fit to dark curve, default 314 is for L=200 nm
        m_eff: effective electron mass, default is that o InAs i.e. 0.023
    Parameters
    ----------
    mu : float
        chemical potential in meV
    kDeltaT_G/D : float
        temperature at steep/gradient said in meV
    V_model : numpy array
        bias voltage vector, in V. must be same number of elements as experimental data
    L : float
        length of barrier in nm
    dE : int/float
        Energy step size, in [meV], for Landau-Buttiker integration
    kT : float
        base temperature in meV, temp in K multiplied by Boltzman in meV/K
    U_top : float
        barrier height in meV
    Ampl : float
        scalling factor between experimental and calculated current, determined from fit to dark curve, default 314 is for L=200 nm
    m_eff : float
        effective electron mass, default is that o InAs i.e. 0.023
        
    Returns
    -------
    I_model:Numpy array
        Array containging calculated current for the corresponding voltage
        '''
    k_B=8.617E-2 #meV /K
    kT = T*k_B  #Temp in meV
    m_e = 9.11E-31 #kg
    
    ###Are units correct?? Double ccheck!!#
    h_bar =  1.05E-34 #Js 4.135E-15 # [eVs]
    V_model = V_model*1E3   #V to mV
    C=2/3*1E-9*np.sqrt(2*m_eff*m_e/(h_bar)**2/(6.24E21)) #conversion constant looping through the potential, see e.g. eq. 7, 1 A is 6-24E18 electrons per second
    
    E = np.arange(mu-10*kT, U_top + 10*kT, dE)# energy integration range, in meV
    I_model = np.zeros(len(V_model))
    # V_model = -V_model
    for i,V in enumerate(V_model):
        f_G=1./(1+np.exp((E-(mu-V))/(kT+kDeltaT_G)))# gradient Fermi fcn. DeltaT is the increase in temp relative base temp, V is in meV
        f_S=1./(1+np.exp((E-mu)/(kT+kDeltaT_S))) # steep Fermi fcn 
        #transmission prob.
        T_barr = np.append(np.exp(-L*C*(U_top-E[E<U_top])**(3/2)/(U_top+V))*(E[E<U_top]), (E[E>=U_top])) #Eq. 7 in Peters manual
        I_model[i] = Ampl*np.sum(T_barr*(f_S-f_G))*dE  #
    return I_model

def Average_IV_sweeps(counter, n_sweeps, I_raw):
    '''
    Parameters
    ----------
    counter : int
        Number of V_h settings.
    n_sweeps : int
        number of total IV sweeps per V_h.
    I_raw : array
        Array of all current sweeps collected.

    Returns
    -------
    I : array
    I_std: array
    '''
    I = np.zeros((I_raw.shape[0],counter))      #array for storing the averaged currents
    I_std = np.zeros((I_raw.shape[0],counter)) 
    for n in range(counter):
        I_temp = I_raw[:,n*n_sweeps:(n+1)*n_sweeps]
        I[:,n] = np.sum(I_temp, axis=1)/n_sweeps
        I_std[:,n] = np.std(I_temp, axis=1)
    return I, I_std

def SumOfSquaredError(model, exp):
    '''
    Calculate SSE between model and experimental data. The dimensions of model and exp must be the same
    Parameters
    ----------
    model : array
        calcualted by model
    exp : array
        experimentally measured

    Returns
    -------
    SSE : array       
    ''' 
    warnings.filterwarnings("ignore") #ignore warnings
    SSE = np.sum((exp - model)**2) #Calculate SE
    return SSE

def fit_WKB_to_data(I_exp, kDeltaT_G, kDeltaT_S, mu, Ampl, L, V_model = np.linspace(-20, 20, 201)):
    '''
 
    Parameters
    ----------
    I_exp : 1D array
        measured current to fit after
    kDeltaT_G : 1D array/float/int
        range of delta temperatures on gradient side to fit for 
    kDeltaT_S : 1D array/float/int
        range of delta temperatures on steep side to fit for 
    mu : float
        Chemical potential.

    Returns
    -------
    I_model : 1D array
        modelled current that minimises error.
    fit_kDT_G : float
        Delta temp on gradient side.
    fit_kDT_S : float
        Delta temp on steepside.
    Error : 1D/2D array
        All errors attempted.

    '''
    # I_model = np.zeros(len(I_exp))
    fit_kDT_G = 0
    fit_kDT_S = 0
    
    if type(kDeltaT_G) and type(kDeltaT_S) == np.ndarray:
        print('fitting both T')
        Error_min = SumOfSquaredError(model = WKB_model(mu, kDeltaT_G[0], kDeltaT_S[0], V_model,Ampl= Ampl, L=L), exp = I_exp) # 1E-10 #initall guess for SumOfSquaredError
        Error = np.zeros((len(kDeltaT_G), len(kDeltaT_S)))
        
        for i, kT_G in enumerate(kDeltaT_G):
            for j, kT_S in enumerate(kDeltaT_S):
                 I_temp = WKB_model(mu = mu, kDeltaT_G = kT_G, kDeltaT_S = kT_S, Ampl=Ampl, L = L, V_model = V_model)
                 Error_temp = SumOfSquaredError(model = I_temp, exp = I_exp)
                 Error[i,j] = Error_temp
                 
                 if Error_temp < Error_min:
                     Error_min = Error_temp
                     # I_model = I_temp
                     fit_kDT_G = kT_G
                     fit_kDT_S = kT_S
                     
    elif type(kDeltaT_S) != np.ndarray:
        print('fitting only gradient side')
        Error_min = SumOfSquaredError(model = WKB_model(mu, kDeltaT_G[0], 0, V_model,Ampl= Ampl, L=L), exp = I_exp) # 1E-10 #initall guess for SumOfSquaredError
        Error = np.zeros(len(kDeltaT_G))
        
        for i, kT_G in enumerate(kDeltaT_G):
            I_temp = WKB_model(mu = mu, kDeltaT_G = kT_G, kDeltaT_S = kDeltaT_S, V_model = V_model, Ampl = Ampl, L=L)
            Error_temp = SumOfSquaredError(model = I_temp, exp = I_exp)
            Error[i] = Error_temp
            
            if Error_temp < Error_min:
                Error_min = Error_temp
                # I_model = I_temp
                fit_kDT_G = kT_G
    
    elif type(kDeltaT_G) != np.ndarray:
        print('fitting only steep side')
        Error_min = SumOfSquaredError(model = WKB_model(mu, 0, kDeltaT_S[0], V_model,Ampl= Ampl, L=L), exp = I_exp) # 1E-10 #initall guess for SumOfSquaredError
        Error = np.zeros(len(kDeltaT_S))
        
        for i, kT_S in enumerate(kDeltaT_S):
            I_temp = WKB_model(mu = mu, kDeltaT_G = kDeltaT_G, kDeltaT_S = kT_S, V_model = V_model, Ampl = Ampl, L=L)
            Error_temp = SumOfSquaredError(model = I_temp, exp = I_exp)
            Error[i] = Error_temp
            
            if Error_temp < Error_min:
                Error_min = Error_temp
                # I_model = I_temp
                fit_kDT_S = kT_S
    else:
        print('make sure correct temperature parameters are provided')  
        
    return [fit_kDT_G, fit_kDT_S, Error]

def fit_A_mu_no_heating(model_type, I_exp, V, A_array, mu_array, L):
    '''
    Parameters
    ----------
    model_type : str
        string defining the type of model to use for calculating transmission, must be one of following:
        'WKB'
        'QM'
    I_exp : 1D array
        measured current to fit after
    V : 1D array
        corresponding voltage bias
    A_array : 1D array
        All values of A to fit for
    mu_array : 1D array
        All values of mu to fit for
    L : float
        length in [nm].

    Returns
    -------
    Error : 2D array
        SSE for all valuese of A and µ checked.

    '''

    
    Error = np.zeros((len(mu_array), len(A_array)))

    if model_type == 'WKB':
        for j, A in enumerate(tqdm(A_array, position=0, leave=True)):        #tqdm is progress bar
            for i, mu in enumerate(mu_array):        
                I_model = WKB_model(mu = mu, kDeltaT_G = 0, kDeltaT_S = 0, V_model = V, Ampl = A, L=L)
                # if Error_temp<Error_min:
                Error_temp = SumOfSquaredError(I_model, I_exp)
                Error[i,j] = Error_temp 

    elif model_type == 'QM':
        print('not ready yet')
    else:
        print('Incorrect model type chosen, use either "WKB" or "QM"')

    return Error


def cut_IV_power_quadrant(I, V):
    '''
    This function shortens an I-V curve so that only the points which lie within the power producing quadrant remains, i.e. 1st or 4th quadrant depending on sign
    ----------
    I : array
        current
    V : array
        voltage

    Returns
    -------
    I_cut : array   
    V_cut : array    
    ''' 
    # I_sc = I[find_nearest(V, 0)]
    V_oc = V[find_nearest(I, 0)]
    if V_oc < 0:        #i.e. 1st quadrant
        start = find_nearest(I, 0)      #index of V_oc
        stop = find_nearest(V, 0)      #index of I_sc

    elif V_oc > 0:        #i.e. 4th quadrant
        start = find_nearest(V, 0)      #index of I_sc
        stop = find_nearest(I, 0)      #index of V_oc
    
    V_cut = V[start:stop]
    I_cut = I[start:stop]
    
    return I_cut, V_cut

def indices_power_quadrant(I, V):
    '''
    This function finds the index in an I-V curve where the power-producing quadrant starts and ends
    ----------
    I : array
        current
    V : array
        voltage

    Returns
    -------
    Start : int 
    stop : int    
    ''' 
    
    V_oc = V[find_nearest(I, 0)]
    if V_oc <= 0:        #i.e. 1st quadrant
        start = find_nearest(I, 0)      #index of V_oc
        stop = find_nearest(V, 0)      #index of I_sc

    elif V_oc > 0:        #i.e. 4th quadrant
        start = find_nearest(V, 0)      #index of I_sc
        stop = find_nearest(I, 0)      #index of V_oc
       
    return start, stop


def Transmission_WKB(E, V, L, U_top, m_eff=0.023):
    '''
    Calculates transmission based on the exponent in the WKB equation, i.e. putting the prefactor C=1
    Parameters
    ----------
    E : int
        energy [eV]
    V : int
        array of bias voltages
    L : int
        length of barrier 
    U_top : int
        Barrier height [eV], default 0.56 eV
    m_eff : int
        effective mass, default is for InAs

    Returns
    -------
    Transmission : int
        The transmission probability, between 0 and 1
    
    ''' 
    m_e = 9.11E-31 #kg
    h_bar = 1.05E-34 #Js
    # h_bar = 6.582*1E-16     #eVs
    h2m = h_bar**2/(2*m_e*m_eff)*6.24E18 # [eVm^2], 6.24E18 converts from J to eV
    #C=2/3*1E-9*np.sqrt(2*m_eff*m_e/(h_bar)**2/(6.24E21))
    #if E<V:
   #     T=0
    if E<U_top:
        T = np.exp(-4/3*L*np.sqrt(1/h2m)*(U_top-E)**(3/2)/(U_top + V))
    elif E>=U_top:
        T=1
    return T

def Transmission_tunnel(E, V, L, U_top, m_eff=0.023):
    '''
    Calculates transmission based on the exponent in the WKB equation, i.e. putting the prefactor C=1
    Parameters
    ----------
    E : int
        energy [eV]
    V : int
        array of bias voltages
    L : int
        length of barrier 
    U_top : int
        Barrier height [eV], default 0.56 eV
    m_eff : int
        effective mass, default is for InAs

    Returns
    -------
    Transmission : int
        The transmission probability, between 0 and 1
    
    ''' 
    from scipy import special # for airy functions
    
    m_e = 9.11E-31 #kg
    h_bar = 1.05E-34 #Js
    # h_bar = 6.582*1E-16     #eVs
    h2m = h_bar**2/(2*m_e*m_eff)*6.24E18 # [eVm^2], 6.24E18 converts from J to eV

    #Various parameters, see Peters description
    alpha = ((U_top+V)/(L*h2m))**(1/3)
    beta = -(E+V)*(L**2/(h2m*(U_top+V)**2))**(1/3)
    q = np.sqrt((E)/h2m)    # if E is set to always be zero at the bottom of band on gradient side then q = np.sqrt((E+V)/h2m) and k = np.sqrt((E)/h2m)
    k = np.sqrt((E-V)/h2m)    
    #The airy functions use to solve differential equation
    AiB, AiBp, BiB, BiBp = special.airy(beta)
    AiaB, AiaBp, BiaB, BiaBp = special.airy(alpha*L+beta)
    
    m_e = 9.11E-31 #kg
    h_bar = 1.05E-34 #Js
    # h_bar = 6.582*1E-16     #eVs
    h2m = h_bar**2/(2*m_e*m_eff)*6.24E18 # [eVm^2], 6.24E18 converts from J to eV
    #C=2/3*1E-9*np.sqrt(2*m_eff*m_e/(h_bar)**2/(6.24E21))
    C = 1/np.pi*(4*q*k*alpha**2)/(alpha**2*np.sqrt(alpha*L+beta)+q**2/np.sqrt(alpha*L+beta))*((alpha*AiBp)**2+(k*AiB)**2)**2
    
    if E<U_top:
        T = C*np.exp(-4/3*L*np.sqrt(1/h2m)*(U_top-E)**(3/2)/(U_top + V))
    elif E>=U_top:
        T=1
    return T

def Transmission_QM(E, V, L, U_top, m_eff=0.023):
    '''
    Parameters
    ----------
    E : int
        energy [eV]
    V : int
        array of bias voltages
    L : int
        length of barrier 
    U_top : int
        Barrier height [eV], default 0.56 eV
    m_eff : int
        effective mass, default is for InAs

    Returns
    -------
    Transmission : int
        The transmission probability, between 0 and 1
    
    ''' 
    from scipy import special # for airy functions

    m_e = 9.11E-31 #kg
    h_bar = 1.05E-34 #Js
    # h_bar = 6.582*1E-16     #eVs
    h2m = h_bar**2/(2*m_e*m_eff)*6.24E18 # [eVm^2], 6.24E18 converts from J to eV

    #Various parameters, see Peters description
    alpha = ((U_top+V)/(L*h2m))**(1/3)
    beta = -(E+V)*(L**2/(h2m*(U_top+V)**2))**(1/3)
    q = np.sqrt((E)/h2m)    # if E is set to always be zero at the bottom of band on gradient side then q = np.sqrt((E+V)/h2m) and k = np.sqrt((E)/h2m)
    k = np.sqrt((E-V)/h2m)    
    #The airy functions use to solve differential equation
    AiB, AiBp, BiB, BiBp = special.airy(beta)
    AiaB, AiaBp, BiaB, BiaBp = special.airy(alpha*L+beta)
    #Calculate transmisison
    G = (1j*alpha*BiaBp+q*BiaB)*(alpha*AiBp+1j*k*AiB)-(1j*alpha*AiaBp+q*AiaB)*(alpha*BiBp+1j*k*BiB)
    T = 4*q*k*alpha**2/(np.pi**2*abs(G)**2)
    if np.isnan(T):
        T= 0
    return T

def fermi_dirac(mu, E, V, kDeltaT_G, kDeltaT_S, T = 77):
    '''
    Parameters
    ----------
    mu : int
        Chemical potential [eV]
    E : array
        energy [eV]
    V : int
        array of bias voltages
    T : int
        ambient temperature [K], default is nitrogen 77K
    Returns
    -------
    f_G: int
    f_S: int
        Fermi-dirac for energy range E
    
    ''' 
    k_B=8.617E-5 #eV /K
    kT = T*k_B  #Temp in eV

    f_G=1./(1+np.exp((E-(mu-V))/(kT+kDeltaT_G)))# gradient Fermi fcn. DeltaT is the increase in temp relative base temp, V is in meV
    f_S=1./(1+np.exp((E-(mu))/(kT+kDeltaT_S)))# steep Fermi fcn. DeltaT is the increase in temp relative base temp
    return f_G, f_S

def Landauer_Buttiker(T, f_S, f_G, A, dE):
    '''
    Calculates current usingn Landauer-Buttiker theory, based on the transmission probability and fermi-dirac distributions
    Parameters
    ----------
    T : array
        Transmission probability
    f_S : array
        fermi-dirac probability distribution for steep (right) side
    f_G : array
        fermi-dirac probability distribution for gradient (left) side
    A : int
        Scalling factor 
    dE : int
        Stepsize for integration
    Returns
    -------
    I_model : Array
        calculated current 
    ''' 
    e = 1.60217E-19 #[C]
    h_eV = 4.1356E-15   #[eV/S]
    I_model = A*np.sum(T*(f_S-f_G))*dE  #(e/h_eV)*
    return I_model

def find_Isc_Voc(I, V):
    '''
    Finds the short circuit current and open circuit voltage for IV data
    Parameters
    ----------
    I : array
        current
    V : array
        voltaqge

    Returns
    -------
    I_sc : float
    V_oc : float
    ''' 
    N = I.shape[1]   #number of IV curves
    I_sc = np.zeros(N)
    V_oc = np.zeros(N)
    for i in range(N):
        I_sc[i] = I[:,i][find_nearest(V, 0)]
        V_oc[i] = V[find_nearest(I[:,i], 0)]
    
    return I_sc, V_oc

def Power(I, V, start = 0, stop = 0):
    '''
    Calculates the power P=IV (not absolute value) for given I-V curve
    ----------
    I : array
        current
    V : array
        voltaqge
    Returns
    -------
    P : array

    ''' 
    if stop != 0:
        V = V[start:stop]
        P = np.zeros((len(V), I.shape[1]))
        for i in range(I.shape[1]):
            P[:,i] = I[:,i][start:stop]*V
    else:
        P = np.zeros((len(V), I.shape[1]))
        for i in range(I.shape[1]):
            P[:,i] = I[:,i]*V
    return P

def Power_cut(I, V):
    '''
    Calculates the power P=IV only in the power producing quadrant 
    ----------
    I : array
        current
    V : array
        voltaqge
    Returns
    -------
    P : array
    V : array
    ''' 
    start, stop = indices_power_quadrant(I, V) 
    I = np.abs(I[start:stop])
    V = np.abs(V[start:stop])
    P = I*V
    return P, V

def fill_factor(I, V):
    '''
    Calculates the Fill factor for given I-V curve
    ----------
    I : array
        current
    V : array
        voltaqge
    Returns
    -------
    FF : Array
    ''' 
    FF = np.zeros(I.shape[1])
#    A = np.zeros(I.shape[1])
#    I_int = np.zeros(I.shape[1])
    for i in range(I.shape[1]):
        V_oc = V[find_nearest(I[:,i], 0)]        #check value of V_oc to know wether power is in 1st or 4th quadrant
        I_sc = I[:,i][find_nearest(V, 0)]
        if V_oc <= 0:        #i.e. 1st quadrant
            start = find_nearest(I[:,i], 0)      #index of V_oc
            stop = find_nearest(V, 0)      #index of I_sc
    
        elif V_oc > 0:        #i.e. 4th quadrant
            start = find_nearest(V, 0)      #index of I_sc
            stop = find_nearest(I[:,i], 0)      #index of V_oc
        else:
            print('wtf')
        #cut out relevant segment
        x = V[start:stop]
        y = I[:,i][start:stop]
            
        I_int = np.trapz(y, x)        
        A = abs(I_sc*V_oc)       
        FF[i] = abs(I_int/A)
        
    return FF

