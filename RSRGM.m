function [psi, theta, theta_star, s] = RSRGM(A, K, lambda, beta, tau, n_iter)
% RSRGM is the main algorithm described in Fig.2.  
%
% Inputs:
%   A: adjacent matrix of PPI network under consideration.
%   K: maximum number of possible functional units.
%   lambda: rate parameter of exponential distribution. The default value
%   is 1.
%   beta: penalization parameter for smooth regularizer. The default value
%   is 1.
%   tau: threshold parameter for obtaining functional units. The default
%   value is 0.3.
%   n_iter: the number of iterations limited in RSRGM. The default value is
%   150.
%
% Outputs:
%   psi: group-protein preference matrix.
%   theta: protein-group membership matrix.
%   theta_star: resultant protein-group indication matrix.
%   s: value of the objective function of (7).


    if nargin < 6
        n_iter = 150;
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
        error('Yor need input adjacent matrix A and maximum possible number K of functional units.');
    end
    
    n = size(A,1);
    D = diag(sum(A,2));
    
    theta = rand(n,K); %Initialize matrix theta randomly.
    psi = rand(K,n); %Initialize matrix psi randomly.
    
    for i  = 1: n_iter
        % Update theta according to Equation (8).
        theta = theta.*(  ((A./(theta*psi+eps))*psi'+ beta*A*psi') ./ (ones(n,n)*psi' +lambda + beta*D*theta +eps) );   
        
        % Update psi according to Equation (9).
        psi = psi.*( (theta'*( A./(theta*psi+eps) ) + beta*theta'*A) ./ (theta'*ones(n,n)+lambda+ beta*psi*D +eps) );
    end
    
    % Calculate the  value s of the objective function (7).  
    s = - sum( sum(A.*log(theta*psi+eps) ) )  + sum(sum(theta*psi)) + lambda*sum(sum(theta)) + lambda*sum(sum(psi)) + 0.5*beta* ( trace(theta'*D*theta) + trace(psi*D*psi') - 2*trace(psi*A*theta) );
    
    % Obtain the resultant protein-group indication matrix theta_star
    % according to Equation (9).
    theta_star = theta;
    theta_star(theta_star >= tau) =1;
    theta_star(theta_star < tau) =0;
    
    % Delete the columns of theta_star that contain at most two elements of
    % 1.
    theta_star (:,sum(theta_star)<= 2) = [];
    
    
    
