%% Direct time integration - Newmark method
% This function simulates the dynamic response of the structure according to
% direct integration, according to a Newmark linear scheme.
% The excitation is a linearly chirped signal at node (4,5,0).
% The response is measured by a three-axial accelerometer at node (4,15,3).

% INPUTS :
%   - time_prop :    Time properties of the simulation              (array)[s]
%   - elements_All : Structure with all initial elements            (struct)[/]
%   - nodes_All :    Exhaustive list of nodes and their DOFs        (array)[/]
%   - M_S :          Structral mass matrix                          (matrix)[/]
%   - C :            Damping matrix                                 (matrix)[/[]
%   - K_S :          Structral stiffness matrix                     (matrix)[/]


% OUTPUTS :
%   - q : Displacement vector according to each DOF, at each time   (array)[/]
%   - q_M_norm : M-norm vector of the dynamic response vector       (array)[/]


function [q,q_M_norm] = newmark(time_prop, elements_All, nodes_All, M_S, C, K_S)


% Initialising the timer
t_newmark = tic;

% Extracting the total number of DOFs
N_DOFs = nodes_All(end,end);

% Extracting the discretisation level
N_elem = numel(fieldnames(elements_All.('Beam1_elements')));

% Extraction of time step and duration
t_max = time_prop(1);           % Duration of the simulation
h = time_prop(2);               % Time step (precision) of the simulation
t = 0:h:t_max;                  % Defining the time counter 

%% Initialisation of external shaker force signal
% TIME BEHAVIOR

% Linearly chirped signal
f_at_0 = 0;             % Initial frequency     (int)[Hz]
t_final = 40;           % Final chirping time   (int)[s]
f_at_t_final = 10;      % Final frequency       (int)[Hz]
phase = 270;            % Phase shift           (int)[°]      

% Constructing the discretised chirped signal
signal = chirp(t,f_at_0,t_final,f_at_t_final,'linear',phase);


% Step signal
% t_step = t_max/4;             % Time of the step variation   (float)[s]
% signal = [zeros(1,fix(t_step/h) + 1) ones(1,fix((t_max - t_step)/h))];

% Dirac signal
% t_dirac = t_max/4;            % Time of the dirac peak        (float)[s]
% signal = [zeros(1,fix(t_dirac/h) + 1) 1 zeros(1,fix((t_max - t_dirac)/h) - 1)];

% Null signal
% signal = zeros(1,fix(t_max/h)+1);

% AMPLITUDE BEHAVIOR
g_force = 250;          % Signal amplitude      (int)[N]

% Defining DOF number corresponding to the excitation
% Excitation is width-oriented at node (4,5,0)
% This node is the final node of the second beam, according to
% the indexation followed during the discretisation

% First, extract the node number, according to 
% final node of last element of beam number 2
shakerNode_nbr = elements_All.('Beam2_elements').(['Element' num2str(N_elem)]).nodeFin_nbr;

% Knowing the number of the node, 
% and as the DOFs are stored in the array nodes_All list
% knowing we have a Y excitation, it corresponds then to the DOF number
shakerDOF_nbr = nodes_All(shakerNode_nbr,6);

% Initialisation of the vector of excitation of constant amplitude 
% (negative because opposed to Y axis)
g = [
    zeros(shakerDOF_nbr-1,1);
    -g_force; 
    zeros(nodes_All(end,end)-shakerDOF_nbr,1)
    ];

% Initialisation of the external force signal array
% It has as rows as the number of DOFs, and as columns as time discretisation
% The element (shakerDOF_nbr,1) will have to be replaced by the value at the
% column corresponding to the time considered
p =  g * signal;


%% Extraction of accelerometer measurements
% Accelerometer is at node (4,15,3). This is the final node of the last element
% of the beam number 16, according to the indexation followed during the 
% discretisation.

% Extracting the corresponding node number
accNode_nbr = elements_All.('Beam16_elements').(['Element' num2str(N_elem)]).nodeFin_nbr;

% Extracting the corresponding 3 translational DOFs from the nodes properties list
accDOFs_nbr = nodes_All(accNode_nbr,5:7);



%% Integration coefficients
gamma = 1/2;
beta = 1/4;


% Initialising the response vectors
% This accounts for null initial conditions as well
q      = zeros(N_DOFs,length(t)); q_st     = zeros(N_DOFs,length(t));
q_dot  = zeros(N_DOFs,length(t)); q_dot_st = zeros(N_DOFs,length(t));
q_ddot = zeros(N_DOFs,length(t));

% Initialisation of the M-norm vector, allowing to study global convergence
q_M_norm = zeros(1,t_max/h);


% Iteration matrix
S = M_S + gamma * h * C + beta * h^2 * K_S;

for t = h:h:t_max
   n = fix(t/h);
   
   % Prediction
   q_dot_st(:,n+1) = q_dot(:,n) + (1 - gamma) * h * q_ddot(:,n);
   q_st(:,n+1) = q(:,n) + h * q_dot(:,n) + (.5 - beta) * h^2 * q_ddot(:,n);
   
   % Computation of accelerations
   p(shakerDOF_nbr,1) = p(shakerDOF_nbr,n + 1);
   q_ddot(:,n+1) = S\(p(:,1) - C * q_dot_st(:,n+1) - K_S * q_st(:,n+1));
   
   % Correction
   q_dot(:,n+1) = q_dot_st(:,n+1) + h * gamma * q_ddot(:,n+1);
   q(:,n+1) = q_st(:,n+1) + h^2 * beta * q_ddot(:,n+1);
   
   % Computing the M-norm
   q_M_norm(fix(t/h)+1) = q(:,fix(t/h)+1)' * M_S * q(:,fix(t/h)+1);
end




tElapsed = toc(t_newmark);
fprintf(['\nDirect integration, Newmark method\n' ... 
'Duration %i s, step %.3f s: %.2f s elapsed\n'],t_max,h,tElapsed);

pl = input('\n Do you want to plot the dynamic response graphs? Y/N\n','s');
if isempty(pl)
    return
elseif strcmp(pl,'y') == 1 || strcmp(pl,'Y') == 1
    figure('Name',...
        ['Direct integration - Newmark method.' ...
        'Time properties : t_max = ' num2str(t_max)...
        ' and h = ' num2str(h)]...
        ,'NumberTitle','off','Color','white'...
        ,'units','normalized','outerposition',[0 0 1 1]);
    title('Dynamic response according to the mode displacement method');
    
    t = 0:h:t_max;
    subplot(2,3,1);
    plot(t,q(accDOFs_nbr(1),fix(t/h)+1));
    title('X-dimension');xlabel('Time [s]'); ylabel('Displacement [m]');
    

    subplot(2,3,2);
    plot(t,q(accDOFs_nbr(2),fix(t/h)+1));
    title('Y-dimension') ;xlabel('Time [s]'); ylabel('Displacement [m]');
    
    subplot(2,3,3);
    plot(t,q(accDOFs_nbr(3),fix(t/h)+1));
    title('Z-dimension'); xlabel('Time [s]'); ylabel('Displacement [m]');
    
    subplot(2,3,4:6);
    plot(t,q_M_norm(fix(t/h)+1));
    title('M-norm'); xlabel('Time [s]');
end



end