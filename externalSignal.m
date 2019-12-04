%% This function computed the external shaker signal
% According to the choice made for the external signal,
% this function outputs a vector whose rows are the external force 
% for each degree of freedom, and the columns are the values of this force
% for each time step.

% INPUTS = 
%   - elements_All  : Structure containing all elements      (struct)[/]
%   - nodes_All     : Array containing all nodes properties (matrix)[/]
%   - time_prop     : Time properties of the simulation     (array)[s]
%   - force_type    : String choosing the type of force     (str)[/]

% OUTPUTS
%   - p             : Vector of external shaker signal      (matrix)[N]
function p = externalSignal(elements_All, nodes_All, time_prop, force_type)

% Extracting the discretisation level
N_elem = numel(fieldnames(elements_All.('Beam1_elements')));
nbr_DOFs = nodes_All(end,end);

% Defining evaluation time step and duration
t_max = time_prop(1);           % Duration of the simulation
h = time_prop(2);               % Time step (precision) of the simulation
t = 0:h:t_max;

%% Initialisation of external shaker force signal
% TIME BEHAVIOR
if strcmp(force_type,'chirp') == 1
    % Linearly chirped signal
    f_at_0 = 0;             % Initial frequency     (int)[Hz]
    t_final = 40;           % Final chirping time   (int)[s]
    f_at_t_final = 10;      % Final frequency       (int)[Hz]
    phase = 270;            % Phase shift           (int)[°]      
    % Constructing the discretised chirped signal
    signal = chirp(t,f_at_0,t_final,f_at_t_final,'linear',phase);

elseif strcmp(force_type,'step') == 1
    % Step signal
    t_step = t_max/4;             % Time of the step variation   (float)[s]
    signal = [zeros(1,fix(t_step/h) + 1) ones(1,fix((t_max - t_step)/h))];

elseif strcmp(force_type,'dirac') == 1
    % Dirac signal
    t_dirac = t_max/4;            % Time of the dirac peak        (float)[s]
    signal = [zeros(1,fix(t_dirac/h) + 1) 1 zeros(1,fix((t_max - t_dirac)/h) - 1)];

elseif strcmp(force_type,'null') == 1
    % Null signal
    signal = zeros(1,fix(t_max/h)+1);
end

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
% knowing we have a X excitation, it corresponds then to the DOF number
shakerDOF_nbr = nodes_All(shakerNode_nbr,5);

% Initialisation of the vector of excitation of constant amplitude 
% (negative because opposed to X axis)
g = [
    zeros(shakerDOF_nbr-1,1);
    -g_force; 
    zeros(nbr_DOFs-shakerDOF_nbr,1)
    ];

% Initialisation of the external force signal array
% It has as rows as the number of DOFs, and as columns as time discretisation
% The element (shakerDOF_nbr,1) will have to be replaced by the value at the
% column corresponding to the time considered
p =  g * signal;










end