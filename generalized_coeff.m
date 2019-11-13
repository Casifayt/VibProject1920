function [mu,gamma,beta] = generalized_coeff(mode_ordre,K_S,M_S,C,V)

x = V(:,mode_ordre);

mu =    x' * M_S * x;
gamma = x' * K_S * x;
beta =  x' *  C  * x;


end