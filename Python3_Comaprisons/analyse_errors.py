# Python36 script to analyse errors
# Author : R. Guasch
# Project : Pre-Doctoral examination, 2D Pincell case study

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm

##-------------------------------------------------------------- Function definitions ---------------------------------------------------------##

#---------------------------------------------------------- Reactivity error/difference analysis -----------------------------------------------#
# Define the relative error function
def error_reactiv(Keffs, Keff_ref):
    """
    Compute the reactivity difference between the reference Keff and the computed Keffs
    """
    errors = []
    for Keff in Keffs:
        errors.append((1/Keff_ref-1/Keff)*1e5)
    return errors

def plot_error_matrix_NXT(errors_matrix, densur, nangle, method):
    fig, ax = plt.subplots()
    cax = ax.matshow(errors_matrix, cmap=cm.bwr)
    fig.colorbar(cax, label='Reactivity difference (pcm)')
    ax.set_title(f'Reactivity difference between NXT {method} and MATLAB {method}', pad=40)
    ax.set_xticks(np.arange(len(nangle)))
    ax.set_yticks(np.arange(len(densur)))
    ax.set_yticklabels(densur)
    ax.set_xticklabels(nangle)
    ax.set_xlabel('Number of angles')
    ax.set_ylabel('Lines density (lines/cm)')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig(f'errors_matrix_NXT_{method}.png')
    plt.show()

def plot_error_matrix_SYBILT(errors_matrix, angular_quadrature_points, number_of_segments):
    fig, ax = plt.subplots()
    cax = ax.matshow(errors_matrix, cmap=cm.coolwarm)
    fig.colorbar(cax, label='Reactivity difference (pcm)')
    ax.set_title('Reactivity difference between SYBILT CP and MATLAB CP', pad=20)
    ax.set_xticks(np.arange(len(angular_quadrature_points)))
    ax.set_yticks(np.arange(len(number_of_segments)))
    ax.set_yticklabels(number_of_segments)
    ax.set_xticklabels(angular_quadrature_points)
    ax.set_xlabel('Number of angular quadrature points')
    ax.set_ylabel('Number of spatial quadrature points')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig('errors_matrix_SYBT_CP.png')
    plt.show()

#---------------------------------------------------------- MOC quadrature order analysis -------------------------------------------------#

def compute_error_order(MOC_MATLAB_KEFF, KEFFS, order_factor):
    errors = []
    for order in MOC_MATLAB_KEFF.keys():
        print(f"order factor {order_factor}")
        print(f"comparing order {order} with {str(int(int(order)*order_factor))}")
        if str(int(int(order)*order_factor)) in KEFFS.keys():
            errors.append((1/MOC_MATLAB_KEFF[order]-1/KEFFS[str(int(int(order)*order_factor))])*1e5)
    return errors

def compute_error_to_MATLAB_nmu4(MOC_MATLAB_KEFF, KEFFS):
    errors = []
    for order in KEFFS.keys():
        errors.append((1/MOC_MATLAB_KEFF-1/KEFFS[order])*1e5)
    return errors


def plot_error_order(orders_list, errors_list, quadrature_methods):
    fig, ax = plt.subplots()
    markers = ['o', 's', 'D']
    colors = ['b', 'r', 'g']
    for i in range(len(orders_list)):
        ax.plot(orders_list[i], errors_list[i], label=f'{quadrature_methods[i]}', marker=markers[i], color=colors[i], markersize=5, linestyle="--")
    ax.set_title(f'Difference in the reactivity between MATLAB and DRAGON5', pad=20)
    ax.set_xlabel('MCCGT angular quadrature order')
    ax.set_ylabel('Reactivity difference (pcm) to MATLAB MOC, nmu=4')
    ax.legend(loc="best")
    ax.grid()
    plt.savefig(f'errors_order_LCMD_vs_CACB_GAUS.png')
    plt.show()

#---------------------------------------------------------- End function definitions -------------------------------------------------------#


#---------------------------------------------------------- Main ---------------------------------------------------------------------------#


densur = [80, 160, 320, 640, 1280] # lines densities
nangle = [10, 20, 40, 80, 160, 320, 640] # number of equidistant angles in [0, pi/2)
NXT_Cp_Keffs = [[1.172364, 1.171430, 1.171299, 1.171542, 1.171460, 1.171481, 1.171451],
                [1.171676, 1.171540, 1.171481, 1.171528, 1.171350, 1.171490, 1.171480],
                [1.171739, 1.171522, 1.171466, 1.171492, 1.171524, 1.171469, 1.171489],
                [1.171930, 1.171672, 1.171636, 1.171622, 1.171628, 1.171606, 1.171597],
                [1.172009, 1.171690, 1.171674, 1.171670, 1.171669, 1.171658, 1.171661]] # lines are densur, columns are nangle 

PIJ_MATLAB_Keff = 1.171670

errors_matrix = [error_reactiv(NXT_Cp_Keffs[i], PIJ_MATLAB_Keff) for i in range(len(densur))]
print(errors_matrix[3][2])
errors_matrix = np.array(errors_matrix)

# plot 2d map of error matrix where x-axis is densur and y-axis is nangle
    
plot_error_matrix_NXT(errors_matrix, densur, nangle, 'CP')


# SYBILT CP

iqua2 = [7, 14, 28, 32]
nseg = [2, 4, 8, 10]
SYBILT_Cp_Keffs = [[1.172514, 1.172548, 1.172554, 1.172567], # rows are nseg, columns are iqua2
                   [1.171676, 1.171620, 1.171622, 1.171627],  
                   [1.171677, 1.171622, 1.171623, 1.171622], 
                   [1.171677, 1.171623, 1.171624, 1.171624]]

errors_matrix_SYB = [error_reactiv(SYBILT_Cp_Keffs[i], PIJ_MATLAB_Keff) for i in range(len(nseg))]
print(errors_matrix_SYB)
plot_error_matrix_SYBILT(errors_matrix_SYB, iqua2, nseg)

# NXT + MCCGT : MOC
MOC_MATLAB_Keff = 1.17161863 
NXT_MOC_nmu4_Keffs =[[1.171366, 1.171391, 1.171420, 1.171423, 1.171421, 1.171422, 1.171422],
                  [1.171396, 1.171422, 1.171443, 1.171450, 1.171451, 1.171451, 1.171452], 
                  [1.171422, 1.171445, 1.171469, 1.171476, 1.171477, 1.171477, 1.171477], 
                  [1.171551, 1.171574, 1.171599, 1.171605, 1.171606, 1.171606, 1.171606], 
                  [1.171600, 1.171623, 1.171647, 1.171653, 1.171654, 1.171655, 1.171655]] 

errors_matrix_NXT_MOC = [error_reactiv(NXT_MOC_nmu4_Keffs[i], MOC_MATLAB_Keff) for i in range(len(densur))]
print(errors_matrix_NXT_MOC[3][1])
print(errors_matrix_NXT_MOC[3][2])
plot_error_matrix_NXT(errors_matrix_NXT_MOC, densur, nangle, 'MOC')


# nmu, quadrature orders
MOC_MATLAB_KEFF = {"2": 1.16461143, "3": 1.17094507, "4": 1.17161863}
LCMD_Keffs = {"2": 1.164567, "3": 1.170901, "4": 1.171574}
CACB_Keffs = {"2": 1.132414, "3": 1.145952, "4": 1.153416, "6": 1.161091, "8": 1.164802, "12": 1.168143, "16": 1.169534, "32": 1.171025}
GAUS_Keffs = {"2": 1.132016, "3": 1.147758, "4": 1.155658, "6": 1.163227, "8": 1.166649, "12": 1.169515, "16": 1.170582, "32": 1.171460}




CACB_error = compute_error_to_MATLAB_nmu4(MOC_MATLAB_KEFF["4"], CACB_Keffs)
GAUS_error = compute_error_to_MATLAB_nmu4(MOC_MATLAB_KEFF["4"], GAUS_Keffs)
LCMD_error = compute_error_to_MATLAB_nmu4(MOC_MATLAB_KEFF["4"], LCMD_Keffs)


LCMD_orders = [2, 3, 4]
CACB_orders = [2, 3, 4, 6, 8, 12, 16, 32]
GAUS_orders = [2, 3, 4, 6, 8, 12, 16, 32]

orders = [LCMD_orders, CACB_orders, GAUS_orders]
errors = [LCMD_error, CACB_error, GAUS_error]
plot_error_order(orders, errors, ["LCMD","CACB","GAUS"])
"""
errors_MOC_LCMD = compute_error_order(MOC_MATLAB_KEFF, LCMD_Keffs, 1)
errors_MOC_CACB = compute_error_order(MOC_MATLAB_KEFF, CACB_Keffs, 1)
errors_MOC_GAUS = compute_error_order(MOC_MATLAB_KEFF, GAUS_Keffs, 1)

errors_MOC_CACB_2 = compute_error_order(MOC_MATLAB_KEFF, CACB_Keffs, 2)
errors_MOC_GAUS_2 = compute_error_order(MOC_MATLAB_KEFF, GAUS_Keffs, 2)

errors_MOC_CACB_4 = compute_error_order(MOC_MATLAB_KEFF, CACB_Keffs, 4)
errors_MOC_GAUS_4 = compute_error_order(MOC_MATLAB_KEFF, GAUS_Keffs, 4)

errors_MOC_CACB_8 = compute_error_order(MOC_MATLAB_KEFF, CACB_Keffs, 8)
errors_MOC_GAUS_8 = compute_error_order(MOC_MATLAB_KEFF, GAUS_Keffs, 8)

print(errors_MOC_GAUS) 
print(errors_MOC_CACB)   
print(errors_MOC_GAUS_2)
print(errors_MOC_CACB_2)
print(errors_MOC_GAUS_4)
print(errors_MOC_CACB_4)

print(errors_MOC_GAUS_8)
print(errors_MOC_CACB_8)
"""