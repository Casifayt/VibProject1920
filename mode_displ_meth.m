%% Mode displacement method
% This function simulates the dynamic response of the structure according to
% modal expansion.
% The excitation is a linearly chirped signal at node (4,5,0).
% The response is measured by a three-axial accelerometer at node (4,15,3).

% INPUTS :
%   - time_prop :    Time properties of the simulation              (array)[s]
%   - elements_All : Structure with all initial elements            (struct)[/]
%   - nodes_All :    Exhaustive list of nodes and their DOFs        (array)[/]
%   - V :            Mode shapes matrix                             (matrix)[/]
%   - M_S :          Structral mass matrix                          (matrix)[/]
%   - wColInV :      Corresponding row for sorted frequencies in V  (array)[/]
%   - w_tr :         Truncated array of sorted frequencies          (array)[rad]
%   - damp_ratios :  Damping ratios of the truncated modes          (array)[/]

% OUTPUTS :
%   - q : Displacement vector according to each DOF, at each time   (array)[/]
%   - q_M_norm : M-norm vector of the dynamic response vector       (array)[/]

function [q,q_M_norm] = mode_displ_meth(time_prop, elements_All, nodes_All,...
    V, M_S, wColInV, w_tr, damp_ratios)

% Initialising the timer
t_displ_meth = tic;

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
% signal = chirp(t,f_at_0,t_final,f_at_t_final,'linear',phase);

% Step signal
% t_step = t_max/4;             % Time of the step variation   (float)[s]
% signal = [zeros(1,fix(t_step/h) + 1) ones(1,fix((t_max - t_step)/h))];

% Dirac signal
% t_dirac = t_max/4;            % Time of the dirac peak        (float)[s]
% signal = [zeros(1,fix(t_dirac/h) + 1) 1 zeros(1,fix((t_max - t_dirac)/h) - 1)];

% Null signal
signal = zeros(1,fix(t_max/h)+1);

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



%% Computation of the simulated response vector

% Extracting the truncation used for the natural frequencies,
% and that will still be used here to truncate the modal expansion
tr = length(w_tr);

% Initialisation of response array
% Column c is the response at time (c-1) * h
% Row r is the response of DOF number r
q = zeros(nodes_All(end,end),t_max/h);

% Initialisation of the M-norm vector, allowing to study global convergence
q_M_norm = zeros(1,t_max/h);

% Initialisation of the inhomogeneous term array
% Row r is the term for mode r
phi = zeros(tr, fix(t_max/h) + 1 );

% Initialisation of the heaviside function vector
heavi = zeros(1, fix(t_max/h) + 1);

for t = 0:h:t_max
    index = fix(t/h)+1;
    
    % Restituting the external force vector as a column vector, at the
    % current time of investigation
    p(shakerDOF_nbr,1) = p(shakerDOF_nbr,index);
        
    for mode = 1:tr
        % Extracting corresponding mode frequency
        omega = w_tr(mode);
        % Extracting corresponding damping ratio
        epsilon = damp_ratios(mode);
        % Defining  damped omega
        omega_damp = omega * sqrt(1 - epsilon^2);
        % Extracting the corresponding mode shape
        mode_shape = V(:,wColInV(mode));
        % Constructing the corresponding generalized mass
        mu = mode_shape' * M_S * mode_shape;
        % Constructing the independant term of the equation
        phi(mode,index) = mode_shape' * p(:,1) / mu;
    end
    
    % Constructing the heaviside function
    if t == 0
        heavi(index) = 0;
    else
        heavi(index) = exp(-epsilon * omega * t) * ...
            sin(omega_damp * t) / omega_damp;
    end
end

prod_conv = zeros(mode,length(phi(mode,:))+length(heavi)-1);

for mode = 1:tr
    prod_conv(mode,:) = conv(phi(mode,:),heavi(:));
end

% Computation
for t = 0:h:t_max
    sumk = 0;                       % Initialisation of the sum
    for mode = 1:tr                 % Sum for each truncated mode
        % Extracting corresponding mode frequency
        omega = w_tr(mode);
        % Extracting corresponding damping ratio
        epsilon = damp_ratios(mode);
        % Defining  damped omega
        omega_damp = omega * sqrt(1 - epsilon^2);
        % Extracting the corresponding mode shape
        mode_shape = V(:,wColInV(mode));
        % Defining the exponential factor
        e = exp(-epsilon * omega * t);
        % Free parameter of the differential equation
        B = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Constructing the solution %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Homogeneous part
        eta_h = e * B * sin(omega_damp * t);
        
        % Inhomogeneous part        
        % Integrating the convolution product
        eta_p = trapz(prod_conv(mode,1:fix(t/h)+1));
        
        % Constructing coefficients
        eta = eta_h + eta_p;
        
        % Summing the series after multiplying by the mode shape
        sumk = sumk + eta * mode_shape;
    end
    
    % Storing the resulting response for the time considered
    % in the dynamic response vector
    q(:,fix(t/h)+1) = sumk;
    % Computing the M-norm
    q_M_norm(fix(t/h)+1) = q(:,fix(t/h)+1)' * M_S * q(:,fix(t/h)+1);
end


tElapsed = toc(t_displ_meth);
fprintf(['\nMode displacement method, %i terms\n' ... 
'Duration %i s, step %.3f s: %.2f s elapsed\n'],tr,t_max,h,tElapsed);

pl = input('\n Do you want to plot the dynamic response graphs? Y/N\n','s');
if isempty(pl)
    return
elseif strcmp(pl,'y') == 1 || strcmp(pl,'Y') == 1
    figure('Name',...
        ['Mode Displacement method truncated to ' num2str(tr) ...
        ' modes. Time properties : t_max = ' num2str(t_max)...
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