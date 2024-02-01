import math
import matplotlib.pyplot as plt
import time

from matplotlib.lines import Line2D
from random import randint, uniform
from typing import List


def compute_T(b):
    return randint(1, 4*b)


def compute_T_hat(T, z):
    eps = (2 * randint(0, 1) - 1) * z
    T_hat = max(T + eps, 1)
    return int(T_hat)


def compute_r0(T_hat, price_list):
    return price_list.index(min(price_list[:T_hat])) + 1 


def compute_r1(price_list, M_p):    
    Pr_r_list = [price_list[i]/(i+1) for i in range(0, M_p)]
    Pr1_r1 = min(Pr_r_list)
    r1 = Pr_r_list.index(Pr1_r1)
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
    
    wr = wi # randint(2,wi)
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
    :param int z: standard deviation
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
    :param Listp[int] price_list: Price list
    :return: ri, Mp and c_opt
    """

    r0 = compute_r0(T_hat, price_list)
    Mp = compute_M_p(price_list)
    r1 = compute_r1(price_list, Mp)
    r2 = compute_r2(price_list, lmbd, r0, r1, Mp)
    r3 = compute_r3(price_list, lmbd, r1, Mp)

    if T_hat >= Mp:
        if r2 > T:
            ri = T
        else:
            ri = price_list[r2-1]
    else:
        if r3 > T:
            ri = T
        else:
            ri = price_list[r3-1]

    ri2 = max(price_list[r2-1] / Mp, price_list[r3-1]/r3)
    comp_ratio = ri / compute_OPT_t(ri, Mp)
    c_opt = price_list[r1-1] / r1

    return comp_ratio, ri2, Mp, c_opt


def plot_multi_lines(
    name,
    data,
    x_label,
    y_label,
    color_list,
    legend_list
):

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


def computeAlg_with_margins_and_plot(Experiments, P, trials, B):
    results_comp_lst = [[] for i in range(Experiments)]
    results_relative_lst = [[] for i in range(Experiments)]
    
    lmbd_list = [i/P for i in range(1, P + 1)]
    lmbd_range = len(lmbd_list)
    print(lmbd_list)
    
    for e in range(Experiments):
        if equilibria :
            price_list = equilibria_generate_price_list(B)
        else:
            price_list = generate_random_price_list(B, z=1)
        print(e)
        print(price_list) 

        competitive_ratio_list = [0] * lmbd_range
        relative_comp_ratio_list = [0] * lmbd_range
        for lmb in range(lmbd_range): 
            lmbd = lmbd_list[lmb]
            competitive_ratio_sum = 0
            relative_comp_ratio_sum = 0

            for i in range(trials):
                T=compute_T(B)
                T_hat = T
                ri, ri2, Mp, c_opt = multiagent_ski_rental(T, T_hat, lmbd, price_list)

                competitive_ratio = ri2
                relative_comp_ratio = (competitive_ratio - 1)/(c_opt - 1)

                competitive_ratio_sum += competitive_ratio
                relative_comp_ratio_sum += relative_comp_ratio

            competitive_ratio_avg = competitive_ratio_sum / trials
            competitive_ratio_list[lmb] = competitive_ratio_avg

            relative_comp_ratio_avg = relative_comp_ratio_sum / trials
            relative_comp_ratio_list[lmb] = relative_comp_ratio_avg
    
            results_comp_lst[e].append((lmb,competitive_ratio_avg))
            results_relative_lst[e].append((lmb,relative_comp_ratio_avg))            
   
    plot_multi_lines(f'competitive_ratio_equilibria_{equilibria}_B_{B}_trials_{trials}_lambda_{P}_eperiments_{Experiments}', \
                     results_comp_lst,
                     x_label="Lambda",
                     y_label="Competitive Ratio",
                     color_list= ['red', 'green', 'blue', 'black', 'orange', 'magenta', 'yellow', 'cyan', 'gray', 'brown'],
                     legend_list= ['price' + str(i + 1) for i in range(Experiments)])

    plot_multi_lines(f'relative_competitive_ratio_equilibria_{equilibria}_B_{B}_trials_{trials}_lambda_{P}_eperiments_{Experiments}', \
                     results_relative_lst,
                     x_label="λ",
                     y_label="Relative Competitive Ratio",
                     color_list= ['red', 'green', 'blue', 'black', 'orange', 'magenta', 'yellow', 'cyan', 'gray', 'brown'],
                     legend_list= ['price' + str(i + 1) for i in range(Experiments)])


def computeAlg_rel_comp_ratio_ws_margins_and_plot(P, trials, B):

    legend_list = ['avg_consistency', 'avg_upperbound']
    results_lst = [[] for i in range(len(legend_list))]
    
    lmbd_list = [i/P for i in range(1, P + 1)]
    lmbd_range = len(lmbd_list)

    if equilibria:
        price_list = equilibria_generate_price_list(B)
    else:
        price_list = generate_random_price_list(B, z=0.1)
    print(price_list) 

    for lmb in range(lmbd_range): 
        lmbd = lmbd_list[lmb]
        ratio_sum = [0 for i in range(len(legend_list))]

        for i in range(trials):
                                
            T = compute_T(B)
            T_hat = T
            ri, ri2, Mp, c_opt = multiagent_ski_rental(T, T_hat, lmbd, price_list)
            ratio_sum[0] += ri
            ratio_sum[1] += ri2

        ratio_avg = [r_sum / trials for r_sum in ratio_sum]
        
        for k in range(len(legend_list)):
            str_ratio_avg = ', '+str(ratio_sum[k])
        print(f"{lmb}{str_ratio_avg}")

        for k in range(len(legend_list)):
            results_lst[k].append((lmb, ratio_avg[k]))

    plot_multi_lines(f'all_ratio_equilibria_{equilibria}_B_{B}_trials_{trials}_lambda_{P}_',
                     results_lst,
                     x_label="Parameter λ",
                     y_label="Competitive Ratio",
                     color_list=['red', 'blue'],
                     legend_list=legend_list)


if __name__ == "__main__":
    equilibria = True
    random = False

    trials = 100
    B = 100
    P = 100

    # Experiments = 10
    # computeAlg_with_margins_and_plot(Experiments, P, trials, B)
    
    computeAlg_rel_comp_ratio_ws_margins_and_plot(P, trials, B)






