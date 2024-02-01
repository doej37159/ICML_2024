import numpy as np
from typing import List

from utils import compute_r0, compute_r1, compute_M_p, compute_r2, compute_r3, generate_equilibria_price_list, \
    generate_random_price_list, compute_T, compute_OPT_t, plot_multi_data


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

    eps = np.random.normal(0, z)
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


def compute_alg(
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

    lmbd_list = [1, 0.2]
    lmbd_range = len(lmbd_list)
    zz_list = [1, 0.9, 0.5]

    for zz in zz_list:
        results_lmbd_lst = [[] for i in range(len(lmbd_list))]
        results_ratio_lst = [[] for i in range(len(lmbd_list))]
        print("z, comp_rat_l_1, comp_rat_l_0.2")
        if equilibria:
            price_list = generate_equilibria_price_list(B)
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

if __name__ == "__main__":
    equilibria = False
    random = False

    P = 250
    trials = 1000
    B = 100

    compute_alg(P, trials, B)
