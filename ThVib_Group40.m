%% MECA0029 - Theory of Vibration
% Analysis of the dynamic behaviour of a truss bridge
% Academic year 2019 - 2020
%
% Casimir FAYT
% ULiege - Aerospatial Engineering

clear all; clc;
t = tic; format long;
tStart = t;

%% General property of whole system
N_dof = 6;              % Number of DOFs of the system  [/]


%% Material properties
% Steel beams

E = 210e9;              % Young's Modulus               [Pa] 
nu = .3;                % Poissons's ratio              [/]
G = .5*E/(1+2*nu);      % Shear modulus                 [Pa]  
rho = 7.8e3;            % Material density              [kg.m^-3] 
% Density for steel extracted from reference book p.390

% All material properties are stored in a structure
mat_prop = struct('E',E,'nu',nu,'rho',rho,'G',G);



%% Initialisation of the beams
addpath Discretisation
print = 0;             % Prints details of function (1) or not (0)
beams_All = beams_initialisation(print);


%% Discretisation of the beams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_elem = 5;            % Number of elements per beam %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print = 0;             % Prints details of function (1) or not (0)
elements_All = discretisation(beams_All,N_elem,print);


%% Construction of the exhaustive nodes list
nodes_All = nodes_list_construction(elements_All,N_dof);


%% Construction of localisation matrix
locel = locel_matrix_init(elements_All, nodes_All);


%% Construction of rotation matrices
rot_matrices = rot_mat_init();


%% Construction of structural matrices
[K_S, M_S] = struct_mat_init...
    (elements_All, rot_matrices, mat_prop, nodes_All, locel);


%% Resolution of eigenvalue system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tr = 8;         % Number of truncated freq (> 1) %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fct = input('\nPlot of the deformed structure? Y/N\n','s');
if strcmp(fct,'y') == 1 || strcmp(fct,'Y') == 1
    plot = 1;
else
    plot = 0;
end

[fEx,f_tr,wEx,w_tr,wColInV,V,D] = eigenSystem(K_S,M_S,tr,plot);



%% Plot of deformed bridge
addpath Plot_functions
if plot == 1
    plotScript(beams_All,elements_All, tr, wColInV, V, locel)
end


%% Construction of damping matrix
[C,damp_ratios] = damp_mat(K_S,M_S,w_tr);


%% Dynamic response
t_max = 40;               % Total duration of response simulation (int)[s]
h = 1e-2;                 % Time step for simulation              (float)[s]
time_prop = [t_max h];

% Possibility to choose the type of excitation signal
forces = {'null','step','dirac','chirp'};
% 1 for no excitation, 2 for a step excitation, 
% 3 for a dirac peak, 4 for the chirp signal
force_type = forces{4};

% Initialisation of the external shaker
p = externalSignal(elements_All, nodes_All, time_prop, force_type);


%% Modal superposition methods
addpath Modal_superposition

% Mode displacement method
% q = mode_displ_meth...
%     (elements_All, nodes_All, V, M_S, wColInV, w_tr, damp_ratios, p, time_prop);


% Mode acceleration method
% q = mode_acc_meth...
%     (elements_All, nodes_All, V, M_S, K_S, wColInV, w_tr, damp_ratios, p, time_prop);


%% Direct time integration - Newmark method
addpath Direct_time_integration
q = newmark...
(time_prop, elements_All, nodes_All, M_S, C, K_S, p);



%% Plot of DFT of displacement vector at accelerometer DOFs
fourier_plot(q, time_prop, elements_All, nodes_All);

%% Model reduction
addpath Model_reduction

% Several choices for the number of degrees of freedom are available
DOF_to_study = [24,48,72];
% Pick one
DOF_r_wanted = DOF_to_study(1);

% Guyan-Irons' method
% [K_red,M_red,C_red, p_red] = guyan_irons...
%       (elements_All, nodes_All, K_S, M_S, C, V, wColInV, tr, p, DOF_r_wanted);

% Craig-Bampton's method
% [K_red,M_red,C_red,p_red] = craig_bampton...
%       (elements_All, nodes_All, K_S, M_S, C, p, DOF_r_wanted);

% eigenSystem_red(K_red,M_red,tr, fEx);

% [q_red,q_M_norm_red] = newmark...
%     (time_prop, elements_All, nodes_All, M_red, C_red, K_red, p_red);



%% Cleaning of the useless variables left
var_2_clean = {'current_beam', 'N_dof', 'nodes_beam_plotting','nodes_elements_ploting',...
    'array_In', 'current_element', 'print','i', 'j','N','beam', 'P1', 'P2', 'pts', 'nu', ...
    'G', 'rho', 'E', 'array_Fin', 'fid', 'eig_freq_cur', 'eig_freq_last', 'tStart', ...
    'tElapsed', 'mode', 'tr', 'f_tr', 'wColInV', 'w_tr','wEx','fEx', 'elements','angle', ...
    'h', 't_max','time_prop','t','f','corr', 'damp_ratios','fct','M_red','K_red',...
    'ans', 'plot', 'forces', 'damp_ratios', 'DOF_R_wanted', 'DOF_to_study',...
    'R_up', 'R_down','total_nbr_dofs','el_nbr','N_dof', 'index', 'var_2_clean'};
clear (var_2_clean{:});