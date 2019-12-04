%% Mode acceleration method
% This function simulates the dynamic response of the structure according to
% modal expansion.
% The excitation is a linearly chirped signal at node (4,5,0).
% The response is measured by a three-axial accelerometer at node (4,15,3).

% INPUTS :
%   - time_prop :    Time properties of the simulation              (array)[s]
%   - elements_All : Structure with all initial elements            (struct)[/]
%   - nodes_All :    Exhaustive list of nodes and their DOFs        (array)[/]
%   - V :            Mode shapes matrix                             (matrix)[/]
%   - M_S :          Structural mass matrix                         (matrix)[/]
%   - K_S :          Structural stiffness matrix                    (matrix)[/]
%   - wColInV :      Corresponding row for sorted frequencies in V  (array)[/]
%   - w_tr :         Truncated array of sorted frequencies          (array)[rad]
%   - damp_ratios :  Damping ratios of the truncated modes          (array)[/]
%   - p :            External excitation force vector               (matrix)[/]

% OUTPUTS :
%   - q : Displacement vector according to each DOF, at each time   (array)[/]

function q = mode_acc_meth(elements_All, nodes_All,...
    V, M_S, K_S, wColInV, w_tr, damp_ratios, p, time_prop)


t_max = time_prop(1);
h = time_prop(2);

% Extracting the truncation used for the natural frequencies,
% and that will still be used here to truncate the modal expansion
tr = length(w_tr);
fprintf('\nMode acceleration method, %i terms\n',tr);

% Initialising the timer
t_acc_meth = tic;

% Extracting the discretisation level
N_elem = numel(fieldnames(elements_All.('Beam1_elements')));


%% Computation of the simulated response vector

% Initialisation of response array
% Column c is the response at time (c-1) * h
% Row r is the response of DOF number r
q = zeros(nodes_All(end,end),t_max/h);

% Initialisation of the M-norm vector, allowing to study global convergence
q_M_norm = zeros(1,t_max/h);

% Initialisation of the inhomogeneous term vector
phi = zeros(tr, fix(t_max/h) + 1 );

% Initialisation of the heaviside function vector
heavi = zeros(tr, fix(t_max/h) + 1);

% Construction of the function suspect to convolution
for t = 0:h:t_max
    index = fix(t/h)+1;
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
        phi(mode,index) = mode_shape' * p(:,index) / mu;
    
        % Constructing the heaviside function
        if t ~= 0
            heavi(mode,index) = exp(-epsilon * omega * t) * ...
                sin(omega_damp * t) / omega_damp;
        end
    end
end

% Initialising the discrete convolution vector
prod_conv = zeros(mode,length(phi(mode,:))+length(heavi)-1);

% Convoluting discrete functions for each mode
for mode = 1:tr
    prod_conv(mode,:) = conv(phi(mode,:),heavi(mode,:)) * h;
end

% Computation
for t = 0:h:t_max
    sum1 = 0;               % Initialisation of the first sum (mode displ.)
    sum2 = 0;               % Initialisation of the second sum (correcting term)
    t_index = fix(t/h) + 1; % Computing time index
        
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
        
        % Parameters of the differential equation according to I.C.
        A = 0;
        B = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Constructing the solution %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Homogeneous part
        eta_h = e * ( A * cos(omega_damp * t) + B * sin(omega_damp * t));
        
        % Inhomogeneous part
        eta_p = prod_conv(mode,t_index);
        
        % Constructing the coefficient
        eta = eta_h + eta_p;
        
        % Summing the series after multiplying by the mode shape
        sum1 = sum1 + eta * mode_shape;
        
        % Second sum
        sum2 = sum2 + phi(mode,t_index) * mode_shape / omega^2;
    end
    
    % Storing the resulting response for the time considered in the response vector
    q(:,t_index) = sum1 + K_S\p(:,t_index) - sum2;
    
    % Computing the M-norm
    q_M_norm(t_index) = q(:,t_index)' * M_S * q(:,t_index);
end

tElapsed = toc(t_acc_meth);
fprintf('Duration %i s, step %.3f s: %.2f s elapsed\n',t_max,h,tElapsed);


pl = input('\nDo you want to plot the dynamic response graphs? Y/N\n','s');
if isempty(pl)
    return
elseif strcmp(pl,'y') == 1 || strcmp(pl,'Y') == 1
    
    %% Extraction of accelerometer measurements DOFs
    % Accelerometer is at node (4,15,3). This is the final node of the last element
    % of the beam number 16, according to the indexation followed during the 
    % discretisation.

    % Extracting the corresponding node number
    accNode_nbr = elements_All.('Beam16_elements').(['Element' num2str(N_elem)]).nodeFin_nbr;

    % Extracting the corresponding 3 translational DOFs from the nodes properties list
    accDOFs_nbr = nodes_All(accNode_nbr,5:7);
    
    
    figure('Name',...
        ['Mode Acceleration method truncated with ' num2str(tr)...
        ' modes. Time properties t_max = ' num2str(t_max)...
        ' and h = ' num2str(h)]...
        ,'NumberTitle','off','Color','white'...
        ,'units','normalized','outerposition',[0 0 1 1]);
    
    title('Dynamic response according to the mode acceleration method');
    
    t = 0:h:t_max;
    
    subplot(2,3,1);
    plot(t,q(accDOFs_nbr(1),fix(t/h)+1));
    title('X-dimension'); xlabel('Time [s]');ylabel('Displacement [m]');
    

    subplot(2,3,2);
    plot(t,q(accDOFs_nbr(2),fix(t/h)+1));
    title('Y-dimension') ; xlabel('Time [s]');ylabel('Displacement [m]');
    
    subplot(2,3,3);
    plot(t,q(accDOFs_nbr(3),fix(t/h)+1));
    title('Z-dimension'); xlabel('Time [s]');ylabel('Displacement [m]');
    
    subplot(2,3,4:6);
    plot(t,q_M_norm(fix(t/h)+1));
    title('M-norm'); xlabel('Time [s]'); ylabel('Amplitude');
end
end