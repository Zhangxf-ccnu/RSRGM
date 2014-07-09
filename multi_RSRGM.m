function [min_psi, min_theta, min_theta_star, min_s] = multi_RSRGM(network, K, lambda, beta, tau, run_times, n_iter)
% multi_RSRGM repeats the entire calculation of RSRGM multiple times and
% chooses the result that gives the lowest value of objective function of
% (7). We also write the cohesive protein complexes revealed by RSRGM to
% file  'cohesive_protein_complexes.txt'  and  the non-cohesive functional
% units revealed by RSRGM to file 'non_cohesive_functional_units.txt'. For
% both two files, each row corresponds an identified functional unit.
%
% Inputs:
%   network: a structure that contains information of PPI network. Filed of
%   adjacent matrix is the adjacent matrix PPI network and filed of protein_list is the
%   name list of correspding proteins.
%   K: maximum number of possible functional units.
%   lambda: rate parameter of exponential distribution. The default value
%   is 1.
%   beta: penalization parameter for smooth regularizer. The default value
%   is 1.
%   tau: threshold parameter for obtaining functional units. The default
%   value is 0.3.
%   run_times: the number of times that we repreat the entire calculation
%   of RSRGM. The default value is 50.
%   n_iter: the number of iterations limited in RSRGM. The default value is
%   150.
%
% Outputs:
%   min_psi: the optimal value of psi that corresponds the result giving
%   the lowest value of objective function of (7).
%   min_theta: the optimal value of theta that corresponds the result giving
%   the lowest value of objective function of (7).
%   min_theta_star: the optimal value of theta_star that corresponds the result giving
%   the lowest value of objective function of (7).
%   min_s: the optimal value of s that corresponds the result giving
%   the lowest value of objective function of (7).


    if nargin < 7
        n_iter = 150;
    end

    if nargin < 6
        run_times = 1;
    end
    
    if nargin < 5
        tau = 0.3;
    end
    
    if nargin < 4
        beta = 1;
    end
    
    if nargin < 3
         lambda = 1; 
    end
    
    if nargin<2
        error('Yor need input the PPI network and maximum possible number K of functional units.');
    end
    
    min_s = inf;
    for i = 1:run_times
        
        disp(['This is the ' num2str(i) '-th run']);
        
        %Main algorithm of RSRGM described in Fig. 2 
        [psi, theta, theta_star, s]=RSRGM(network.adjacent_matrix, K, lambda, beta, tau, n_iter); 
        
        % Choose the result that gives the lowest value of objective
        % function of  (7).
        if s < min_s
            min_psi = psi;
            min_theta = theta;
            min_theta_star = theta_star;
            min_s = s;
        end
  
    end
    
    
    
   
   %%
   % Calculate the  density of funcaional units revealed by RSRGM.
    density = [];
    for i = 1:length(min_theta_star(1,:))
        t = find(min_theta_star(:,i));
        density(i) = sum(sum(network.adjacent_matrix(t,t)))/(length(t)*(length(t)-1));
    end
    
    % The identified functional units with density above 0.1 are viewed as
    % cohesive complexes, and functional units with density below 0.1 are
    % viewed as non-cohesive functional unis.
    low_density_indices = find(density < 0.1);
    high_density_indices = find(density >= 0.1);
    
    % Write the results of cohesive protein complexes to file
    % 'cohesive_protein_complexes.txt', where each row corresponds to an identified  complex. 
    fid = fopen('cohesive_protein_complexes.txt','w');
    for i = 1:length(high_density_indices)
        member_indices = find(min_theta_star(:,high_density_indices(i)));
        for j = 1:length(member_indices)
            fprintf(fid, '%s\t', cell2mat(network.protein_list(member_indices(j))) );
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    % Write the results of non-cohesive functional units to file
    % 'non_cohesive_functional_units.txt', where each row corresponds to an identified non-cohesive functional unit. 
    fid = fopen(['non_cohesive_functional_units.txt'],'w');
    for i = 1:length(low_density_indices)
        member_indices = find(min_theta_star(:,low_density_indices(i)));
        for j = 1:length(member_indices)
            fprintf(fid, '%s\t', cell2mat(network.protein_list(member_indices(j))) );
        end
        fprintf(fid, '\n');
    end
    fclose(fid);