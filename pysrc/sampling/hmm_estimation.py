import numpy as np
from hmmlearn import hmm
from scipy.linalg import expm, logm

from pysrc.services.data_service import load_cattle_prices


def estimate_price_model(dist=[0.5, 0.5], num_iter=5, var="uncon"):
    for _ in range(num_iter):
        s_low, s_high = dist[0], dist[1]

        # Estimate HMM transitions and states
        aic, ll, bic, price_states, sigmas, P = _estimate(s_low, s_high, var)

        # Compute stationary distribution
        dist, sta_price = _stationary_distribution(P, price_states)

        # Annualize transition matrix
        M = _annualize_transition_matrix(P)

    return price_states, M


def _estimate(s_low=0.5, s_high=0.5, var="uncon"):
    # Read data
    price = load_cattle_prices()
    log_price = np.log(price)

    # Estimating the model
    Q = log_price.reshape(log_price.shape[0], 1)

    np.random.seed(123)

    best_score = best_model = None
    n_fits = 500

    for idx in range(n_fits):
        if var == "uncon":
            model = hmm.GaussianHMM(
                n_components=2, random_state=idx, init_params="tmc", params="stmc"
            )
        else:
            model = hmm.GaussianHMM(
                n_components=2,
                random_state=idx,
                init_params="tmc",
                params="stmc",
                covariance_type="tied",
            )

        model.startprob_ = np.array([s_low, s_high])
        model.fit(Q)
        score = model.score(Q)

        if best_score is None or score > best_score:
            best_model = model
            best_score = score

    aic = best_model.aic(Q)
    bic = best_model.bic(Q)
    mus = np.exp(np.ravel(best_model.means_))
    sigmas = np.ravel(np.sqrt([np.diag(c) for c in best_model.covars_]))
    P = best_model.transmat_

    sorted_indices = np.argsort(mus)
    mus_sorted = mus[sorted_indices]
    sigmas_sorted = sigmas[sorted_indices]
    P_sorted = P[sorted_indices][:, sorted_indices]

    if var == "uncon":
        ll = (best_model.aic(Q) - 12) / (-2)
    else:
        ll = (best_model.aic(Q) - 10) / (-2)

    return (aic, ll, bic, mus_sorted, sigmas_sorted, P_sorted)


def _stationary_distribution(transition_matrix, mus):
    eigenvals, eigenvects = np.linalg.eig(transition_matrix.T)
    close_to_1_idx = np.isclose(eigenvals, 1)
    target_eigenvect = eigenvects[:, close_to_1_idx]
    target_eigenvect = target_eigenvect[:, 0]
    stationary_distrib = target_eigenvect / sum(target_eigenvect)
    stationary_price = mus[0] * stationary_distrib[0] + mus[1] * stationary_distrib[1]
    return (stationary_distrib, stationary_price)


def _annualize_transition_matrix(transition_matrix):
    dτ = 1 / 12  # time step
    P = transition_matrix  # probability transition matrix
    M = logm(P) / dτ  # instantenous generator
    P = expm(M)  # Updated Probability transition matrix
    return P
