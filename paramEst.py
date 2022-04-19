# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 12:06:34 2022

@author: mwu
"""

import matplotlib.pyplot as plt
from scipy import optimize, integrate, interpolate
import numpy as np
import pandas as pd
from thresholding import dapi_and_ki67_analysis

def ODEfunc(t,y,params):
    k_G0G1, k_G1G0, k_G1S, k_SG2M, k_G2MG1 = params 
    G0, G1, S, G2M = y    
    react1 = k_G0G1*G0  
    react2 = k_G1G0*G1  
    react3 = k_G1S*G1   
    react4 = k_SG2M*S    
    react5 = k_G2MG1*G2M   
    k_apop = react5/100
    dG0 = react2 - react1 - k_apop*G0                     # [1/hr] dG0/dt 
    dG1 = react1 + 2*react5 - react2 - react3 - k_apop*G1   # [1/hr] dG1/dt 
    dS = react3 - react4 - k_apop*S                      # [1/hr] dS/dt
    dG2M = react4 - react5 - k_apop*G2M                    # [1/hr] dG2/dt
                    
    return [dG0, dG1, dS, dG2M]

# load experimental data
xdata1 = [48]
y_dmso, y_nrg, y_nrg_p38, y_tt10 = dapi_and_ki67_analysis()

# initial conditions
y0 = y_dmso
tspan = [0, 50]
paramsFixed = [1, 1] # fixed reaction rate constants for k_G1G0 and k_G2MG1
p0 = [1/8, 1, 1] # rudimentary "guesses" for reaction rates for k_G0G1, k_G1S, and k_SG2M

bar_label = np.array(["k_G0G1", "k_G1S", "k_SG2M"])

# run parameter estimation
def objectiveFunction(paramsEst,paramsFixed,tspan,xdata,ydata):
    k_G1G0, k_G2MG1 = paramsFixed
    k_G0G1, k_G1S, k_SG2M = paramsEst
    params = [k_G0G1, k_G1G0, k_G1S, k_SG2M, k_G2MG1]
    sol = integrate.solve_ivp(ODEfunc,tspan,y0,args=(params,))
    finterp = interpolate.interp1d(sol.t,sol.y[i,:].T, fill_value="extrapolate")
    error = ydata - finterp(xdata) 
    return error

height_list = []
std_list = []

# loop through to find the parameter estimates
for i in range(0, 3):
    # dmso
    result_dmso = optimize.least_squares(objectiveFunction,p0,args=(paramsFixed,tspan,xdata1,y_dmso[i]))
    paramsEst_dmso = result_dmso.x
    
    # nrg
    result_nrg = optimize.least_squares(objectiveFunction,p0,args=(paramsFixed,tspan,xdata1,y_nrg[i]))
    paramsEst_nrg = result_nrg.x
    
    # nrg + p38
    result_nrg_p38 = optimize.least_squares(objectiveFunction,p0,args=(paramsFixed,tspan,xdata1,y_nrg_p38[i]))
    paramsEst_nrg_p38 = result_nrg_p38.x
        
    # tt10
    result_tt10 = optimize.least_squares(objectiveFunction,p0,args=(paramsFixed,tspan,xdata1,y_tt10[i]))
    paramsEst_tt10 = result_tt10.x

    # covariance dmso
    covariance_dmso = result_dmso.cost*np.linalg.inv((np.matmul([[result_dmso.jac[0][i].T]],[[result_dmso.jac[0][i]]])))/len(y_dmso)
    n = 10;
    paramEnsemble_dmso = np.random.multivariate_normal([paramsEst_dmso[i]],covariance_dmso, n)
    paramsStd_dmso = np.std(paramEnsemble_dmso,axis=0)
    
    # covariance nrg
    covariance_nrg = result_nrg.cost*np.linalg.inv((np.matmul([[result_nrg.jac[0][i].T]],[[result_nrg.jac[0][i]]])))/len(y_nrg)
    n = 10;
    paramEnsemble_nrg = np.random.multivariate_normal([paramsEst_nrg[i]],covariance_nrg,n)
    paramsStd_nrg = np.std(paramEnsemble_nrg,axis=0)
    
    # covariance nrg+p38
    covariance_nrg_p38 = result_nrg_p38.cost*np.linalg.inv((np.matmul([[result_nrg_p38.jac[0][i].T]],[[result_nrg_p38.jac[0][i]]])))/len(y_nrg_p38)
    n = 10;
    paramEnsemble_nrg_p38 = np.random.multivariate_normal([paramsEst_nrg_p38[i]],covariance_nrg_p38,n)
    paramsStd_nrg_p38 = np.std(paramEnsemble_nrg_p38,axis=0)
    
    # covariance tt10
    covariance_tt10 = result_tt10.cost*np.linalg.inv((np.matmul([[result_tt10.jac[0][i].T]],[[result_tt10.jac[0][i]]])))/len(y_tt10)
    n = 10;
    paramEnsemble_tt10 = np.random.multivariate_normal([paramsEst_tt10[i]],covariance_tt10,n)
    paramsStd_tt10 = np.std(paramEnsemble_tt10,axis=0)
    
    # creating arrays for both the parameter estimates and standard deviation
    height = np.array([paramsEst_dmso[i], paramsEst_nrg[i], paramsEst_nrg_p38[i], paramsEst_tt10[i]])
    height_list.append(height)
    
    std = np.array([paramsStd_dmso[0], paramsStd_nrg[0], paramsStd_nrg_p38[0], paramsStd_tt10[0]])
    std_list.append(std)
    
# converting the list to an numpy array
height_list = np.asarray(height_list)
std_list = np.asarray(std_list)
    
# converting the numpy array to the pandas dataframe
height_df = pd.DataFrame(height_list)
std_df = pd.DataFrame(std_list)

height_df.loc[3] = [1, 1, 1, 1]
std_df.loc[3] = [0, 0, 0, 0]

# plotting parameter estimation bar graph
treatments = ["k_G0G1", "k_G1S", "k_SG2M", "k_G2MG1"]
x_axis = np.arange(len(treatments))
w = 0.1

plt.bar(x_axis -w*1.5, height_df.loc[:, 0], width=w, label = 'dmso', color = '#d7191c')
plt.bar(x_axis -w/2, height_df.loc[:, 1], width=w, label = 'nrg', color = '#fdae61')
plt.bar(x_axis +w/2, height_df.loc[:, 2], width=w, label = 'nrg + p38', color = '#CF9FFF')
plt.bar(x_axis +w*1.5, height_df.loc[:, 3], width=w, label = 'tt10', color = '#2c7bb6')

plt.errorbar(x_axis -w*1.5, height_df.loc[:, 0], std_df.loc[:, 0], ls = 'none', color = 'black')
plt.errorbar(x_axis -w*1.5, height_df.loc[:, 1], std_df.loc[:, 1], ls = 'none', color = 'black')
plt.errorbar(x_axis -w*1.5, height_df.loc[:, 2], std_df.loc[:, 2], ls = 'none', color = 'black')
plt.errorbar(x_axis -w*1.5, height_df.loc[:, 3], std_df.loc[:, 3], ls = 'none', color = 'black')

plt.xticks(x_axis, treatments)
plt.legend()
plt.title("Estimated Parameters Bar Plot", fontweight = "bold")
plt.ylabel("Estimated Parameters")
plt.ylim([0, 1.2])
plt.show()

## plotting CM cell cycle concetration vs time graphs with parameter estimates
# parameter estimates for each experimental condition
p0_dmso = height_df.loc[0:2, 0]
p0_nrg = height_df.loc[0:2, 1]
p0_nrg_p38 = height_df.loc[0:2, 2]
p0_tt10 = height_df.loc[0:2, 3]

curve_label = np.array(["G0", "G1", "S", "G2/M"])

# looping through to plot
for i in range(0, 4):
    # run parameter estimation
    def objectiveFunction(paramsEst,paramsFixed,tspan,xdata,ydata):
        k_G1G0, k_G1S = paramsFixed
        k_G0G1, k_SG2M, k_G2MG1 = paramsEst
        params = [k_G0G1, k_G1G0, k_G1S, k_SG2M, k_G2MG1]
        sol = integrate.solve_ivp(ODEfunc,tspan,y0,args=(params,))
        finterp = interpolate.interp1d(sol.t,sol.y[i,:].T, fill_value="extrapolate")
        error = ydata - finterp(xdata) 
        return error
    
    # dmso
    result_dmso = optimize.least_squares(objectiveFunction,p0_dmso,args=(paramsFixed,tspan,xdata1,y_dmso[i]))
    paramsEst_dmso = result_dmso.x
    
    # nrg
    result_nrg = optimize.least_squares(objectiveFunction,p0_nrg,args=(paramsFixed,tspan,xdata1,y_nrg[i]))
    paramsEst_nrg = result_nrg.x
    
    # nrg + p38
    result_nrg_p38 = optimize.least_squares(objectiveFunction,p0_nrg_p38,args=(paramsFixed,tspan,xdata1,y_nrg_p38[i]))
    paramsEst_nrg_p38 = result_nrg_p38.x
        
    # tt10
    result_tt10 = optimize.least_squares(objectiveFunction,p0_tt10,args=(paramsFixed,tspan,xdata1,y_tt10[i]))
    paramsEst_tt10 = result_tt10.x

    # plotting parameter estimates and curves
    params = [paramsEst_dmso[0], paramsFixed[0], paramsFixed[1], paramsEst_dmso[1], paramsEst_dmso[2]]
    sol_dmso = integrate.solve_ivp(ODEfunc,tspan,y0,args=(params,))
    plt.plot(xdata1,y_dmso[i],'o', markerfacecolor = '#d7191c', color = '#d7191c')
    l1, = plt.plot(sol_dmso.t,sol_dmso.y[i,:].T,'-', color = '#d7191c')
    
    params = [paramsEst_nrg[0], paramsFixed[0], paramsFixed[1], paramsEst_nrg[1], paramsEst_nrg[2]]
    sol_nrg = integrate.solve_ivp(ODEfunc,tspan,y0,args=(params,))
    plt.plot(xdata1, y_nrg[i], 'o', markerfacecolor = '#fdae61', color = '#fdae61')
    l2, = plt.plot(sol_nrg.t,sol_nrg.y[i,:].T,'-', color = '#fdae61')
    
    params = [paramsEst_nrg_p38[0], paramsFixed[0], paramsFixed[1], paramsEst_nrg_p38[1], paramsEst_nrg_p38[2]]
    sol_nrg_p38 = integrate.solve_ivp(ODEfunc, tspan, y0, args=(params,))
    plt.plot(xdata1, y_nrg_p38[i], 'o', markerfacecolor = '#CF9FFF', color = '#CF9FFF')
    l3, = plt.plot(sol_nrg_p38.t, sol_nrg_p38.y[i, :].T, '-', color = '#CF9FFF')
    
    params = [paramsEst_tt10[0], paramsFixed[0], paramsFixed[1], paramsEst_tt10[1], paramsEst_tt10[2]]
    sol_tt10 = integrate.solve_ivp(ODEfunc, tspan, y0, args=(params,))
    plt.plot(xdata1, y_tt10[i], 'o', markerfacecolor = '#2c7bb6', color = '#2c7bb6')
    l4, = plt.plot(sol_tt10.t, sol_tt10.y[i, :].T, '-', color = '#2c7bb6')
    
    plt.legend([l1, l2, l3, l4], ["dmso", "nrg", "nrg + p38", "tt10"])
    plt.title(curve_label[i], fontweight = "bold")
    plt.xlabel("Time [hr]")
    plt.ylabel("Percent of cells in " + curve_label[i])
    plt.ylim([0, 100])
    plt.show()
