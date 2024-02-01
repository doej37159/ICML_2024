from random import randint
from typing import List

from utils import compute_r0, compute_r1, compute_M_p, compute_r2, compute_r3, generate_equilibria_price_list, \
    generate_random_price_list, compute_T, compute_OPT_t, plot_multi_lines


def compute_T_hat(
    T,
    z
):
    """
    Computes T hat.

    :param T: Number of skying days
    :param z: Standard deviation
    :return: T hat
    """

    eps = (2 * randint(0, 1) - 1) * z
    T_hat = max(T + eps, 1)
    return int(T_hat)


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


def compute_alg_with_margins(
    Experiments,
    P,
    trials,
    B
):
    """

    :param Experiments: Number of experiments
    :param P: Quantization parameter
    :param trials: Number of trials per experiment
    :param B: Maximum price
    :return: void
    """

    results_comp_lst = [[] for i in range(Experiments)]
    results_relative_lst = [[] for i in range(Experiments)]
    
    lmbd_list = [i/P for i in range(1, P + 1)]
    lmbd_range = len(lmbd_list)
    print(lmbd_list)
    
    for e in range(Experiments):
        if equilibria :
            price_list = generate_equilibria_price_list(B)
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


def compute_alg_rel_comp_ratio_ws_margins(
    P,
    trials,
    B
):
    """

    :param P: Quantization parameter
    :param trials: Number of trials
    :param B: Maximum price
    :return: void
    """

    legend_list = ['avg_consistency', 'avg_upperbound']
    results_lst = [[] for i in range(len(legend_list))]
    
    lmbd_list = [i/P for i in range(1, P + 1)]
    lmbd_range = len(lmbd_list)

    if equilibria:
        price_list = generate_equilibria_price_list(B)
    else:
        price_list = generate_random_price_list(B, z=0.1)

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
    
    compute_alg_rel_comp_ratio_ws_margins(P, trials, B)
