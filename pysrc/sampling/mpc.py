import numpy as np
import pandas as pd
import os
from ..services.data_service import load_site_data
from pysrc.sampling import baseline

def mc_samples_constrained(location):
    
    np.random.seed(123)
    # Simulate samples for model unconstrained variance
    num_simulations = 200

    states = {'low': 32.44, 'high': 42.78}
    initial_state = 'high' 

    probability_matrix = {
        'low': {'low': 0.766, 'high': 0.234},
        'high': {'low': 0.046, 'high': 0.954}
    }

    # Number of observations
    num_observations = 200


    for i in range(1, num_simulations + 1):
        
        # Generating the Markov chain for each simulation
        current_state = initial_state
        markov_chain = [states[current_state]]

        for _ in range(num_observations - 1):
            next_state = np.random.choice(list(probability_matrix[current_state].keys()), 
                                        p=list(probability_matrix[current_state].values()))
            markov_chain.append(states[next_state])
            current_state = next_state

        transformed_markov_chain = [1 if price == 32.44 else 2 for price in markov_chain]

        # Creating an index from 1 to 200
        index = range(1, num_observations + 1)

        # Combine index and transformed Markov chain into a DataFrame
        markov_chain_df = pd.DataFrame({'Index': index, 'scenario': transformed_markov_chain})

        # Specify the filename for each CSV file (e.g., mc_1.csv, mc_2.csv, ..., mc_100.csv)
        csv_filename = f'/mc_{i}.csv'
        output=location
        # Output to CSV file
        markov_chain_df.to_csv(output+csv_filename, index=False)
    
    
    return("mc sampling is done")


def mc_samples_unconstrained(location):
    
    np.random.seed(123)

    num_simulations = 200

    states = {'low': 35.7, 'high': 44.3}
    initial_state = 'high' 

    probability_matrix = {
        'low': {'low': 0.706, 'high': 0.294},
        'high': {'low': 0.171, 'high': 0.829}
    }

    # Number of observations
    num_observations = 200


    for i in range(1, num_simulations + 1):
        # Generating the Markov chain for each simulation
        current_state = initial_state
        markov_chain = [states[current_state]]

        for _ in range(num_observations - 1):
            next_state = np.random.choice(list(probability_matrix[current_state].keys()), 
                                        p=list(probability_matrix[current_state].values()))
            markov_chain.append(states[next_state])
            current_state = next_state

        transformed_markov_chain = [1 if price == 35.7 else 2 for price in markov_chain]

        # Creating an index from 1 to 200
        index = range(1, num_observations + 1)

        # Combine index and transformed Markov chain into a DataFrame
        markov_chain_df = pd.DataFrame({'Index': index, 'scenario': transformed_markov_chain})

        # Specify the filename for each CSV file (e.g., mc_1.csv, mc_2.csv, ..., mc_100.csv)
        csv_filename = f'/mc_{i}.csv'
        output=location
        # Output to CSV file
        markov_chain_df.to_csv(output+csv_filename, index=False)

    return("mc sampling is done")



def gdx_files(location,num_sites=78):
    # Load site data
    
    
    (
        zbar_2017,
        z_2017,
        forest_area_2017,
        _,
        _,
        _,
        _,
    ) = load_site_data(num_sites)
        
    
    baseline_fit = baseline.sample(
    model_name="full_model",
    num_sites=num_sites,
    iter_sampling=10**4,
    chains=5,
    seed=1,
    )

    theta = baseline_fit.stan_variable("theta").mean(axis=0)
    gamma = baseline_fit.stan_variable("gamma").mean(axis=0)
    
    
    scale=1e9
    # Computing carbon absorbed in start period
    x0_vals = gamma * forest_area_2017

    x0_vals = x0_vals * scale
    site_z_vals = z_2017 * scale
    zbar_2017 = zbar_2017 * scale

    x0data = pd.DataFrame(x0_vals)
    x0data.columns=['x_2017_78Sites']
    x0data.index = x0data.index + 1
    saveto = os.path.join(location, "X0Data.csv")
    x0data.to_csv(saveto)

    z0data = pd.DataFrame(site_z_vals)
    z0data.columns=['z_2017_78Sites']
    z0data.index = z0data.index + 1
    saveto = os.path.join(location, "Z0Data.csv")
    z0data.to_csv(saveto)

    zbardata = pd.DataFrame(zbar_2017)
    zbardata.columns=['zbar_2017_78Sites']
    zbardata.index = zbardata.index + 1
    saveto = os.path.join(location, "ZbarData.csv")
    zbardata.to_csv(saveto)

    gammadata = pd.DataFrame(gamma)
    gammadata.columns = ['gamma_78Sites']
    gammadata.index = gammadata.index + 1
    saveto = os.path.join(location, "GammaData.csv")
    gammadata.to_csv(saveto)

    thetadata = pd.DataFrame(theta)
    thetadata.columns=['theta_78Sites']
    thetadata.index = thetadata.index + 1
    saveto = os.path.join(location, "ThetaData.csv")
    thetadata.to_csv(saveto)

    return("gdx file is done")