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
%   - K_S :          Structural stiffness matrix                    (matrix)[/]
%   - p :            External excitation force vector               (matrix)[N]


% OUTPUTS :
%   - q : Displacement vector according to each DOF, at each time   (array)[/]
%   - q_M_norm : M-norm vector of the dynamic response vector       (array)[/]


function q = newmark(time_prop, elements_All, nodes_All, M, C,...
    K, p)

fprintf('\nDirect integration, Newmark method\n'); 

% Initialising the timer
t_newmark = tic;

% Extracting the total number of DOFs
N_DOFs = length(K);

% Extracting the discretisation level
N_elem = numel(fieldnames(elements_All.('Beam1_elements')));

% Extraction of time step and duration
t_max = time_prop(1);           % Duration of the simulation
h = time_prop(2);               % Time step (precision) of the simulation
n_time_step = fix(t_max/h);
t = 0:h:t_max;                  % Defining the time counter 


%% Integration coefficients
gamma = .5;
beta = .25;


% Initialisation of the displacement and velocity response vectors
% This accounts for null initial conditions as well
q      = zeros(N_DOFs,length(t)); q_st     = zeros(N_DOFs,length(t));
q_dot  = zeros(N_DOFs,length(t)); q_dot_st = zeros(N_DOFs,length(t));
q_ddot = zeros(N_DOFs,length(t)); 

% Initial acceleration
q_ddot(:,1) = M\p(:,1);

% Initialisation of the M-norm vector, allowing to study global convergence
q_M_norm = zeros(1,t_max/h);

% Iteration matrix
S = M + h * gamma * C + h^2 * beta * K;
inv_S = inv(S);

for n = 2:n_time_step
   
   % Prediction
   q_dot_st(:,n) = q_dot(:,n-1) + (1 - gamma) * h * q_ddot(:,n-1);
   q_st(:,n) = q(:,n-1) + h * q_dot(:,n-1) + (.5 - beta) * h^2 * q_ddot(:,n-1);
   
   % Computation of accelerations
   q_ddot(:,n) = inv_S * (p(:,n) - C * q_dot_st(:,n) - K * q_st(:,n));
   
   % Correction
   q_dot(:,n) = q_dot_st(:,n) + h * gamma * q_ddot(:,n);
   q(:,n) = q_st(:,n) + h^2 * beta * q_ddot(:,n);
   
   % Computing the M-norm
   q_M_norm(n) = q(:,n)' * M * q(:,n);
end

%% Extraction of accelerometer measurements
% Accelerometer is at node (4,15,3). This is the final node of the last element
% of the beam number 16, according to the indexation followed during the 
% discretisation.

% Extracting the corresponding node number
accNode_nbr = elements_All.('Beam16_elements')...
    .(['Element' num2str(N_elem)]).nodeFin_nbr;

% Extracting the corresponding 3 translational DOFs from the nodes properties list
accDOFs_nbr = nodes_All(accNode_nbr,5:7);


tElapsed = toc(t_newmark);
fprintf('Duration %i s, step %.3f s: %.2f s elapsed\n',t_max,h,tElapsed);

pl = input('\nDo you want to plot the dynamic response graphs? Y/N\n','s');
if isempty(pl)
    return
elseif strcmp(pl,'y') == 1 || strcmp(pl,'Y') == 1
    figure('Name',...
        ['Direct integration - Newmark method.' ...
        'Time properties : t_max = ' num2str(t_max)...
        ' and h = ' num2str(h)]...
        ,'NumberTitle','off','Color','white'...
        ,'units','normalized','outerposition',[0 0 1 1]);
    title('Dynamic response according to the Newmark integration scheme');
    
    t = h:h:t_max;
    
    subplot(3,1,1);
    plot(t,q(accDOFs_nbr(1),fix(t/h)));
    title('X-dimension');xlabel('Time [s]'); ylabel('Displacement [m]');

% Below is the code used when comparing response vector from different methods
    
%     for i = 1:length(qMD(accDOFs_nbr(1),:))
%         if q(accDOFs_nbr(2),i) ~= 0 && q(accDOFs_nbr(1),i) > 1e-12
%             eps(i) = abs(q(accDOFs_nbr(1),i) - qMD(accDOFs_nbr(1),i))/q(accDOFs_nbr(1),i);
%         end
%     end
%     eps_avg_X = sum(eps)/length(eps);
%     plot(t,qMD(accDOFs_nbr(1),fix(t/h)) - qMA(accDOFs_nbr(1),fix(t/h)));
%     legend('Newmark - Mode Displacement');%, 'Mode Displacement - Mode Acceleration');
    

    subplot(3,1,2);
    plot(t,q(accDOFs_nbr(2),fix(t/h)));
    title('Y-dimension') ;xlabel('Time [s]'); ylabel('Displacement [m]');

% Below is the code used when comparing response vector from different methods

%     for i = 1:length(qMD(accDOFs_nbr(2),:))
%         if q(accDOFs_nbr(2),i) ~= 0 && q(accDOFs_nbr(2),i) > 1e-12
%             eps(i) = abs(q(accDOFs_nbr(2),i) - qMD(accDOFs_nbr(2),i))/q(accDOFs_nbr(2),i);
%         end
%     end
%     eps_avg_Y = sum(eps)/length(eps);
%     
%     plot(t,qMD(accDOFs_nbr(2),fix(t/h)) - qMA(accDOFs_nbr(2),fix(t/h)));
%     legend('Newmark - Mode Displacement');%, 'Mode Displacement - Mode Acceleration');
    
    subplot(3,1,3);
    plot(t,abs(q(accDOFs_nbr(3),fix(t/h))));
    title('Z-dimension'); xlabel('Time [s]'); ylabel('Displacement [m]');

% Below is the code used for comparing response vector from different methods

%     for i = 1:length(q(accDOFs_nbr(3),:))
%         if q(accDOFs_nbr(3),i) ~= 0 && q(accDOFs_nbr(3),i) > 1e-12
%             eps(i) = abs(q(accDOFs_nbr(3),i) - qMD(accDOFs_nbr(3),i))/q(accDOFs_nbr(3),i);
%         end
%     end
%     eps_avg_Z = sum(eps)/length(eps);
%     fprintf('eps_avg = %e\n',(eps_avg_X + eps_avg_Y + eps_avg_Z)/3);
%     plot(t,qMD(accDOFs_nbr(3),fix(t/h)) - qMA(accDOFs_nbr(3),fix(t/h)));
%     legend('Newmark - Mode Displacement');, 'Mode Displacement - Mode Acceleration');
    

% Below is the code used to compute the M-norm
%     subplot(2,3,4:6);
%     plot(t,q_M_norm(fix(t/h)));
%     title('M-norm'); xlabel('Time [s]');
end
end