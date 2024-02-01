import math
import time

from random import randint, uniform
from typing import List
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


###############################
##### Algorithm Functions #####
###############################


def calculate_wr(
    r, 
    b
):
    """
    Computes wr.
    
    :param r: Number of days
    :param b: Max price
    :return: wr
    """
    
    wr = 1

    while (wr - 1) * min(b, r - 1 + wr) <= r * (b - 1):
        wr += 1

    return wr


def generate_equilibria_price_list(
    b: int
) -> List[int]:
    """
    Generates an equilibria case price list.
    
    :param int b: Max price
    :return: Equilibrium price list
    """

    price_list: List[int] = []
    r = b
    wr = calculate_wr(r, b)

    for j in range(0, 4*b):
        if j == r:
            Pj = j + wr
        else:
            Pj = j + b
        price_list.append(Pj)

    return price_list


def generate_random_price_list(
    b: int,
    z: float
) -> List[int]:
    """
    Generates a random case price list.
    
    :param int b: Max price
    :param int z: standard deviation
    :return: Random price list
    """

    price_list: List[int] = []
    psi = int(uniform(0, z*b))
    for i in range(0, 4*b):
        Pi = i + b - psi
        price_list.append(Pi)
    return price_list


def compute_M_p(price_list):
    """
    Computes the minimum price in the price list.
    
    :param price_list: Price list
    :return: Minimum price
    """
    
    return min(price_list)


def compute_T(b):
    """
    Computes the number of skying days.
    
    :param b: Maximum price
    :return: Number of skying days
    """
    
    return randint(1, 4*b)


def compute_OPT_t(
    t, 
    M_p
):
    """
    Computes the optimum value for a given day t.
    
    :param t: Day t
    :param M_p: Minimum price
    :return: Optimum value for day t
    """
    
    if t < M_p:
        return t
    else:
        return M_p


def compute_r0(
    T_hat, 
    price_list
):
    """
    Computes r0.
    
    :param T_hat: Predicted ski days
    :param price_list: Price list
    :return: r0
    """
    
    return price_list.index(min(price_list[:T_hat])) + 1


def compute_r1(
    price_list, 
    M_p
):
    """
    Computes r1.
    
    :param price_list: Price list
    :param M_p: Minimum price
    :return: r1
    """
    
    Pr_r_list = [price_list[i]/(i+1) for i in range(0, M_p)]
    Pr1_r1 = min(Pr_r_list)
    r1 = Pr_r_list.index(Pr1_r1)

    return r1 + 1


def compute_r2(
    price_list, 
    lmbd, 
    r0, 
    r1, 
    M_p
):
    """
    Computes r2
    
    :param price_list: Price list
    :param lmbd: Lambda
    :param r0: r0
    :param r1: r1
    :param M_p: Minimum price
    :return: r2
    """
    
    r2_start = math.ceil((1 - lmbd) * (r0 - 1) + lmbd * r1) - 1
    Pr2_OPTr2_list = [price_list[i] - compute_OPT_t(i+1, M_p) for i in range(r2_start,len(price_list))]
    diff_Pr2_OPTr2_min = min(Pr2_OPTr2_list)
    r2 = r2_start + Pr2_OPTr2_list.index(diff_Pr2_OPTr2_min)
    Pr2 = price_list[r2]

    if Pr2 > price_list[r1 - 1]:
        print(f"changing {r2} for {r1}")
        r2 = r1

    return r2 + 1


def compute_r3(
    price_list, 
    lmbd, 
    r1, 
    M_p
):
    """
    Computes r3

    :param price_list: Price list
    :param lmbd: Lambda
    :param r1: r1
    :param M_p: Minimum price
    :return: r3
    """
    
    Pr1 = price_list[r1 - 1]
    c_opt = Pr1/r1
    Pr3_r3_min = math.inf
    r3_min = -1

    for p in range(0, len(price_list)):
        r3 = p + 1
        opt_r3 = compute_OPT_t(r3, M_p)
        Pr3 = price_list[p]

        if (Pr3 / opt_r3) <= 1 + ((1 / lmbd) * (c_opt - 1)) and \
            Pr3 / r3 < Pr3_r3_min:
            Pr3_r3_min = Pr3 / r3
            r3_min = r3

    return r3_min


###############################
##### Plot Functions #####
###############################


def plot_multi_data(
    name, 
    data, 
    x_label='Epsilon', 
    y_label='Competitive Ratio'
):
    """

    :param name: Name of the plot
    :param data: Plot data
    :param x_label: Label for x axis
    :param y_label: Label for y axis
    :return: void
    """
    
    fig = plt.figure(figsize=(7, 7))

    # Plot Title
    plt.title(f"Competitive Ratio of Equilibrium", pad=38)

    # Plot labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Plot Data
    x_val_1 = [x[0] for x in data[0]]
    y_val_1 = [x[1] for x in data[0]]
    plt.plot(x_val_1, y_val_1, color='blue')

    x_val_2 = [x[0] for x in data[1]]
    y_val_2 = [x[1] for x in data[1]]
    plt.plot(x_val_2, y_val_2, color='orange')

    plt.grid()

    # Plot Limits
    ax = plt.gca()
    ax.set_ylim([min(min(y_val_1), min(y_val_2)) - 0.2,
                 max(max(y_val_1), max(y_val_2)) + 0.2])

    custom = [Line2D([], [], marker='.', markersize=10, color='blue', linestyle='None'),
              Line2D([], [], marker='.', markersize=10, color='orange', linestyle='None')]

    plt.legend(custom, ['λ = 1', 'λ = 0.20'], bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
               borderaxespad=0, ncol=2, fontsize=10, frameon=False)

    # Save Plot
    fig.savefig(f"outputs/{name + time.strftime('%Y_%m_%d_%H_%M_%S')}")


def plot_multi_lines(
    name,
    data,
    x_label,
    y_label,
    color_list,
    legend_list
):
    """

    :param name: Name of the plot
    :param data: Plot data
    :param x_label: Label for x axis
    :param y_label: Label for y axis
    :param color_list:  Plot color list
    :param legend_list: Plot legend list
    :return: void
    """
    
    fig = plt.figure(figsize=(7, 7))

    # Plot Title
    plt.title(f"Consistency of Equilibrium", pad=38)

    # Plot labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Plot Data
    ymin = []
    ymax = []
    custom = []
    for line in range(len(data)):
        x_val = [x[0] for x in data[line]]
        y_val = [x[1] for x in data[line]]
        ymin.append(min(y_val))
        ymax.append(max(y_val))
        plt.plot(x_val, y_val, color=color_list[line])
        custom.append(Line2D([], [], marker='.', markersize=10, color=color_list[line], linestyle='None'))

    y_lim_min = min(ymin)
    y_lim_max = max(ymax)

    plt.grid()

    # Plot Limits
    ax = plt.gca()
    ax.set_ylim([y_lim_min - 0.3, y_lim_max + 0.3])

    plt.legend(custom, legend_list, bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
               borderaxespad=0, ncol=5, fontsize=10, frameon=False)

    fig.savefig(f"outputs/{name + time.strftime('%Y_%m_%d_%H_%M_%S')}")
