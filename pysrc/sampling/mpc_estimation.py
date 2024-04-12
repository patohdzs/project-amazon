import numpy as np 
import pandas as pd
import os
from hmmlearn import hmm
from scipy.linalg import logm, expm
import warnings

def est(s_low=0.5,s_high=0.5,var='uncon'):
    
    warnings.filterwarnings("ignore", message="KMeans is known to have a memory leak on Windows with MKL, when there are less chunks than available threads. You can avoid it by setting the environment variable OMP_NUM_THREADS=2.")

    
    # need to input initial state prob s_low and s_high
    data_folder = os.getcwd()+"/data/hmc/"
    
    # read data
    df = pd.read_csv(data_folder+"seriesPriceCattle_prepared.csv")
    price = df['price_real_mon_cattle'].values.astype(float) 
    logprice=np.log(price)
    
    
    # estimating the model
    
    Q=logprice.reshape(logprice.shape[0],1)

    np.random.seed(123)

    best_score = best_model = None
    n_fits = 500

    for idx in range(n_fits):

        if var=='uncon':
            model = hmm.GaussianHMM(n_components=2,random_state=idx,init_params='tmc',params='stmc')
        else:
            model = hmm.GaussianHMM(n_components=2,random_state=idx,init_params='tmc',params='stmc',covariance_type='tied')
            
        model.startprob_= np.array([s_low, s_high])
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
    
    if var=='uncon':
        ll = (best_model.aic(Q)-12)/(-2)
    else:
        ll = (best_model.aic(Q)-10)/(-2)
        
    
    return(aic,ll,bic,mus_sorted,sigmas_sorted,P_sorted)




def stationary(transition_matrix,mus):
    
    eigenvals, eigenvects = np.linalg.eig(transition_matrix.T)
    close_to_1_idx = np.isclose(eigenvals,1)
    target_eigenvect = eigenvects[:,close_to_1_idx]
    target_eigenvect = target_eigenvect[:,0]
    stationary_distrib = target_eigenvect / sum(target_eigenvect)
    stationary_price=mus[0]*stationary_distrib[0]+mus[1]*stationary_distrib[1] 
    return(stationary_distrib,stationary_price)


def annual_transition(transition_matrix):
    dτ = 1/12  # time step
    P = transition_matrix  #probability transition matrix
    M = logm(P)/dτ  # instantenous generator
    P = expm(M)  # Updated Probability transition matrix
    return P


def iteration_est(initial_prob, num_iterations=5,var='uncon'):
    
    
    for i in range(num_iterations):
        if i == 0:
            s_low, s_high = initial_prob[0], initial_prob[1]
        else:
            s_low, s_high = sta_dist[0], sta_dist[1]
        
        aic, ll, bic, mus, sigmas, P = est(s_low, s_high, var)
        sta_dist, sta_price = stationary(P, mus)
        annual_P=annual_transition(P)
        
    return aic,ll,bic,mus,sigmas,P,sta_dist,sta_price,annual_P
    
