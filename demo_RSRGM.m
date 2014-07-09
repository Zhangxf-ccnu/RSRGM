% Choose a network from the following four networks to test RSRGM.

database='Gavin';
% database='Krogan'
% database='Collins';
% database='BioGRID';

switch database
    
    case 'Gavin'
        load ./data/Gavin_network.mat
        K = 500;
        lambda = 1;
        beta = 2;
        tau = 0.3;
        run_times = 50;
        n_iter = 150;
        [psi, theta, theta_star, s] = multi_RSRGM(Gavin_network, K, lambda, beta, tau, run_times, n_iter);
        
    case 'Krogan'
        load ./data/Krogan_network.mat
        K = 500;
        lambda = 2;
        beta = 1;
        tau = 0.3;
        run_times = 50;
        n_iter = 150;
        [psi, theta, theta_star, s] = multi_RSRGM(Krogan_network, K, lambda, beta, tau, run_times, n_iter);
        
    case 'Collins'
        load ./data/Collins_network.mat
        K = 500;
        lambda = 2;
        beta = 1;
        tau = 0.3;
        run_times = 50;
        n_iter = 150;
        [psi, theta, theta_star, s] = multi_RSRGM(Collins_network, K, lambda, beta, tau, run_times, n_iter);
         
    case 'BioGRID'
        load ./data/BioGRID_network.mat
        K = 1000;
        lambda =1;
        beta = 1;
        tau = 0.3;
        run_times = 50;
        n_iter = 150;
        [psi, theta, theta_star, s] = multi_RSRGM(BioGRID_network, K, lambda, beta, tau, run_times, n_iter);
          
end


     