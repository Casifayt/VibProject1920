%% This function constructs the damping matrix
% Under proportional assumption giving
%            C = a K + b M
% Where a and b are determined thanks to a condition on the two first
% damping ratios such that they are equal to 1%

% INPUTS
%   - K_S :         Structural stiffness matrix                     (matrix)[/]
%   - M_S :         Structural mass matrix                          (matrix)[/]
%   - w_tr :        List of truncated modes                         (array)[rad]

% OUTPUTS
%   - C :           Damping matrix                                  (matrix)[/]
%   - damp_ratios : List of the damping ratios up to mode nbr tr    (array)[/]
function [C,damp_ratios] = damp_mat(K_S,M_S,w_tr)

% Initialising the coefficients of the system
A = [
    .5*w_tr(1), .5/w_tr(1);
    .5*w_tr(2), .5/w_tr(2);
    ];

% Initialising the independant term
B = [ .01; .01];

% Solving the problem
sol = A\B;

% Constructing the damping matrix
C = sol(1) * K_S + sol(2) * M_S;

% Initialising the array of damping ratios
damp_ratios = zeros(1,length(w_tr));

%fprintf('Damping ceofficients : %e s and %e Hz\n',sol(1),sol(2));
fprintf('\nDamping ratios under proportional assumption : \n');
for i = 1:length(w_tr)
   ratio = .5*(sol(1)*w_tr(i)+sol(2)/w_tr(i)); 
   damp_ratios(i) = ratio;
   fprintf('Mode %i : \x03B5 = %.2f %%\n',i,100*ratio);
end
end