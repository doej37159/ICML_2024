import math
import matplotlib.pyplot as plt
import time
import numpy as np

from matplotlib.lines import Line2D
from random import randint, uniform
from typing import List


def compute_T(b):
    return randint(1, 4*b)


def compute_T_hat(T, z):
    # z is 0.1, 0.2, ... 1 percent
    eps = np.random.normal(0, z)
    T_hat = max(T + eps, 1)
    return int(T_hat)


def compute_r0(T_hat, price_list):
    return price_list.index(min(price_list[:T_hat])) + 1 


def compute_r1(price_list, M_p):    
    Pr_r_list = [price_list[i]/(i+1) for i in range(0, M_p)]
    Pr1_r1 = min(Pr_r_list)
    r1 = Pr_r_list.index(Pr1_r1)
    # r1 is index and r1+1 are days
    return r1 + 1


def compute_M_p(price_list):
    return min(price_list)


def compute_OPT_t(t, M_p):
    if t < M_p:
        return t
    else:
        return M_p


def compute_r2(price_list, lmbd, r0, r1, M_p):
    r2_start = math.ceil((1 - lmbd) * (r0 - 1) + lmbd * r1) - 1
    Pr2_OPTr2_list = [price_list[i] - compute_OPT_t(i+1, M_p) for i in range(r2_start,len(price_list))]
    diff_Pr2_OPTr2_min = min(Pr2_OPTr2_list)
    r2 = r2_start + Pr2_OPTr2_list.index(diff_Pr2_OPTr2_min)
    Pr2 = price_list[r2]
    if Pr2 > price_list[r1 - 1]:
        print(f"changing {r2} for {r1}")
        r2 = r1
    return r2 + 1


def compute_r3(price_list, lmbd, r1, M_p):
    Pr1 = price_list[r1 - 1]
    c_opt = Pr1/r1
    Pr3_r3_min = math.inf
    r3_min = -1
    for p in range(0, len(price_list)):
        # r3 is the current day
        r3 = p + 1
        opt_r3 = compute_OPT_t(r3, M_p)
        Pr3 = price_list[p]
        if (Pr3 / opt_r3) <= 1 + ((1 / lmbd) * (c_opt - 1)) and \
            Pr3 / r3 < Pr3_r3_min:
            # r3 will always be the first day for which the minimum ratio is achieved by using the < condition
            Pr3_r3_min = Pr3 / r3
            r3_min = r3
    return r3_min


def calculate_wr(r, b):
    wi = 1
    while (wi - 1) * min(b, r - 1 + wi) <= r*(b-1):
        wi += 1
    
    wr = wi
    return wr


def equilibria_generate_price_list(
    b: int
) -> List[int]:
    """
    :param int b: Max price
    :return: Price list
    """

    price_list: List[int] = []
    r = b
    wr = calculate_wr(r, b)

    print(f"r = {r+1}, wr = {wr}")

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
    :param int b: Max price
    :param int z: Standard deviation
    :return: Price list
    """

    price_list: List[int] = []
    psi = int(uniform(0, z*b))
    for i in range(0, 4*b):
        Pi = i + b - psi
        price_list.append(Pi)
    return price_list


def multiagent_ski_rental(
    T: int,
    T_hat: int,
    lmbd: float,
    price_list: List[int]
):
    """
    :param int T: Total Ski Days
    :param int T_hat: Predicted Ski Days, T_hat
    :param float lmbd: Lambda
    :param List[int] price_list: Price list
    :return: ri, Mp and c_opt
    """

    r0 = compute_r0(T_hat, price_list)
    Mp = compute_M_p(price_list)
    r1 = compute_r1(price_list, Mp)
    r2 = compute_r2(price_list, lmbd, r0, r1, Mp)
    r3 = compute_r3(price_list, lmbd, r1, Mp)

    if T_hat >= Mp:
        if r2 > T:
            ri = -1
        else:
            ri = r2
    else:
        if r3 > T:
            ri = -1
        else:
            ri = r3
    c_opt = price_list[r1-1] / r1

    return ri, Mp, c_opt


def plot_multi_data(name, data, x_label='Epsilon', y_label='Competitive Ratio'):

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
    fig.savefig(f"outputs/{name+time.strftime('%Y_%m_%d_%H_%M_%S')}")


if __name__ == "__main__":
    equilibria = False
    random = False

    P = 250
    trials = 1000
    B = 100
    lmbd_list = [1, 0.2]
    lmbd_range = len(lmbd_list)
    zz_list = [1, 0.9, 0.5]
    Experiments = 1

    for zz in zz_list:
        results_lmbd_lst = [[] for i in range(len(lmbd_list))]
        results_ratio_lst = [[] for i in range(len(lmbd_list))]
        print("z, comp_rat_l_1, comp_rat_l_0.2")
        if equilibria:
            price_list = equilibria_generate_price_list(B)
        else:
            price_list = generate_random_price_list(B, zz)
        print(price_list)

        for z in range(1, 250):
            competitive_ratio_list = [0] * lmbd_range
            relative_comp_ratio_list = [0] * lmbd_range

            for lmb in range(lmbd_range): 
                lmbd = lmbd_list[lmb]
                competitive_ratio_sum = 0
                competitive_ratio_max = 0
                relative_comp_ratio_sum = 0
                relative_comp_ratio_max = 0
                center = compute_r1(price_list, compute_M_p(price_list))

                for i in range(trials):
                    T = compute_T(B)
                    T_hat = compute_T_hat(T, z-1)
                    ri, Mp, c_opt = multiagent_ski_rental(T, T_hat, lmbd, price_list)
                    OPT = compute_OPT_t(T, Mp)

                    if ri < 0:
                        CA = T
                    else:
                        CA = price_list[ri - 1]

                    competitive_ratio = CA/OPT
                    relative_comp_ratio = (competitive_ratio - 1)/(c_opt - 1)

                    competitive_ratio_sum += competitive_ratio
                    relative_comp_ratio_sum += relative_comp_ratio

                    competitive_ratio_avg = competitive_ratio_sum / trials
                    competitive_ratio_list[lmb] = competitive_ratio_avg

                    relative_comp_ratio_avg = relative_comp_ratio_sum / trials
                    relative_comp_ratio_list[lmb] = relative_comp_ratio_avg
            for i in range(lmbd_range):
                results_lmbd_lst[i].append((z, competitive_ratio_list[i]))

            for i in range(lmbd_range):
                results_ratio_lst[i].append((z, relative_comp_ratio_list[i]))

            print(f"{z},{competitive_ratio_list[0]}, {competitive_ratio_list[1]}, ...")

        plot_multi_data(f'competitive_ratio_equilibria_{equilibria}_B_{B}_trials_{trials}_qu_{P}_z_0_{int(zz*10)}_',\
                        results_lmbd_lst,
                        x_label="Standard Deviation σ", 
                        y_label="Competitive Ratio")
