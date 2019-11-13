%% MECA0029 - Theory of Vibration
% Analysis of the dynamic behaviour of a truss bridge
% Academic year 2019 - 2020
%
% Casimir FAYT
% ULiege - Aerospatial Engineering

clear all; clc; f = fopen('Log.txt','w'); fclose(f);
t = tic; format long;
tStart = t;
logpr('Program launched');

%% General property of whole system
N_dof = 6;              % Number of DOFs of the system  [/]


%% Material properties
% Steel beams

E = 210e9;              % Young's Modulus               [Pa] 
nu = .3;                % Poissons's ratio              [/]
G = .5*E/(1+2*nu);      % Shear modulus                 [Pa]  
rho = 7.8e3;            % Material density              [kg.m^-3] 
% Density for steel extracted from reference book p.390

mat_prop = struct('E',E,'nu',nu,'rho',rho,'G',G);
% All material properties are stored in a structure


%% Initialisation of the beams
logpr('Initialisation of the beams');
print = 0;             % Prints details of function (1) or not (0)
beams_All = beams_initialisation(print);
t = logprt('Beams initialised',t);


%% Discretisation of the beams
N_elem = 2;            % Number of elements per beam
logpr(['Discretisation into ' num2str(N_elem) ' elements']);
print = 0;             % Prints details of function (1) or not (0)
elements_All = discretisation(beams_All,N_elem,print);
t = logprt('Beams discretised',t);


%% Construction of the exhaustive nodes list
logpr('Construction of the exhaustive nodes properties list');
nodes_All = nodes_list_construction(elements_All,N_dof);
t = logprt('Nodes list constructed',t);


%% Construction of localisation matrix
logpr('Construction of localisation matrix');
locel = locel_matrix_init(elements_All, nodes_All);
t = logprt('Localisation matrix constructed',t);


%% Construction of rotation matrices
logpr('Construction of rotation matrices');
rot_matrices = rot_mat_init();
t = logprt('Rotation matrices constructed',t);


%% Construction of structural matrices
logpr('Construction of structural matrices');
[K_S, M_S] = struct_mat_init(elements_All, rot_matrices, mat_prop, nodes_All, locel);
t = logprt('Structural matrices constructed',t);

%% Resolution of eigenvalue system
tr = 8;         % Number of truncated freq (> 1)
logpr(['Resolution of the eigenvalue system,' ...
    'frequencies truncated to ' num2str(tr)]);
[fEx,f_tr,wEx,w_tr,wColInV,V,D] = eigenSystem(K_S,M_S,tr);
t = logprt('Eigenvalue system solved',t);


%% Construction of deformed elements
mode = 5;                  % Sorted number of the eigenmode desired for plot
logpr(['Construction of the deformed elements,'...
    'according to mode number ' num2str(mode)]);
elements_Deformed = deformation(mode,wColInV,V,locel,elements_All);
t = logprt('Deformed elements constructed.',t);
tElapsed = toc(tStart);
logpr(['\nTotal elapsed time so far : ' num2str(tElapsed) 's\n']);


%% Plot of deformed bridge
plotPerso(beams_All, elements_Deformed, tr);


%% Construction of damping matrix
[C,damp_ratios] = damp_mat(K_S,M_S,w_tr);


%% Dynamic response
t_max = 50;
h = t_max*1e-3;
time_prop = [t_max h];

% Mode displacement method
% [q,q_M_norm] = mode_displ_meth(time_prop,elements_All,nodes_All,...
% V, M_S, wColInV,w_tr,damp_ratios);

% Mode acceleration method
% [q,q_M_norm] = mode_acc_meth(time_prop, elements_All,nodes_All,...
%     V, M_S, K_S, wColInV,w_tr,damp_ratios);

% Direct time-integration Newmark method
[q,q_M_norm] = newmark(time_prop, elements_All, nodes_All, M_S, C, K_S);


%% Cleaning of the useless variables left
var_2_clean = {'current_beam', 'N_dof', 'nodes_beam_plotting','nodes_elements_ploting',...
    'array_In', 'current_element', 'print','i', 'j','N','beam', 'P1', 'P2', 'pts', 'nu', ...
    'G', 'rho', 'E', 'array_Fin', 'fid', 'eig_freq_cur', 'eig_freq_last', 'tStart', ...
    'tElapsed', 'mode', 'tr', 'f_tr', 'wColInV', 'w_tr','wEx','fEx', 'elements','angle', ...
    'h', 't_max','time_prop','t','f','corr', 'damp_ratios'...
    'R_up', 'R_down','total_nbr_dofs','el_nbr','N_dof', 'index', 'var_2_clean'};
clear (var_2_clean{:});



%% Functions for the log file
function t =  logprt(str,timer)
% Function to write in log file with timer

% Opening log file
log = fopen('Log.txt', 'a');
if log == -1
  error('Cannot open log file.');
end

tElapsed = toc(timer);

% Printing in log file
fprintf(log, '%s: %s in %f s\n', datestr(now, 'dd/mm/yy'), str, tElapsed);

% Closing the log file
fclose(log);

% Starting the new timer
t = tic;
end

function logpr(str)
% Function to write in log file without timer


% Opening log file
log = fopen('Log.txt', 'a');
if log == -1
  error('Cannot open log file.');
end

% Printing in log file
fprintf(log, '%s: %s\n', datestr(now,'dd/mm/yy'), str);

% Closing the log file
fclose(log);
end