%% Function of model redution using the dynamic Craig Bampton's method
% This function reduces the problem size into a certain number of retained
% nodes. This reductions costs in accuracy but yields a high decrease in
% computation time.

% INPUTS
%   - elements_All : Structure with all elements properties     (structure)[/]
%   - nodes_All :    Array with all nodes and their DOFs        (structure)[/]
%   - K_S :          Structural stiffness matrix                (matrix)[/]
%   - M_S :          Structural mass matrix                     (matrix)[/]
%   - C :            Damping matrix                             (matrix)[/]
%   - p :            External excitation force vector           (matrix)[N]
%   - DOF_r_wanted : Number of degrees of freedom retained      (int)[/]

% OUTPUTS
%   - K_red :       Reduced stiffness matrix                    (matrix)[/]
%   - M_red :       Reduced mass matrix                         (matrix)[/]
%   - C_red :       Reduced damping matrix                      (matrix)[/]
%   - p_red :       Reduced external excitation force vector    (matrix)[N]

function [K_red,M_red,C_red, p_red] = craig_bampton(elements_All, nodes_All,...
    K, M, C, p, DOF_r_wanted)

% Initialising timer
t_CB = tic;

N_elem = numel(fieldnames(elements_All.('Beam1_elements')));
nbr_DOFs = nodes_All(end,end);

% Extraction of remaining nodes
% Remaining nodes are the nodes (0,15,0); (4,15,0);(0,15,3); (4,15,3)
% The two first are extreme nodes of beam number 4 in current indexation
% The two last are extreme nodes of beam number 16 in current indexation


% Extracting the number of the nodes
nR_nodes_nbr(1:4) = [
elements_All.('Beam4_elements').('Element1').nodeIn_nbr, ...
elements_All.('Beam4_elements').(['Element' num2str(N_elem)]).nodeFin_nbr,...
elements_All.('Beam16_elements').('Element1').nodeIn_nbr,...
elements_All.('Beam16_elements').(['Element' num2str(N_elem)]).nodeFin_nbr...
];

if DOF_r_wanted ~= 24
    nR_nodes_nbr(5:8) = [
    elements_All.('Beam3_elements').('Element1').nodeIn_nbr,...
    elements_All.('Beam3_elements').(['Element' num2str(N_elem)]).nodeFin_nbr...
    elements_All.('Beam15_elements').('Element1').nodeIn_nbr,...
    elements_All.('Beam15_elements').(['Element' num2str(N_elem)]).nodeFin_nbr...
    ];
end

if DOF_r_wanted == 72
    nR_nodes_nbr(9:12) = [
    elements_All.('Beam2_elements').('Element1').nodeIn_nbr,...
    elements_All.('Beam2_elements').(['Element' num2str(N_elem)]).nodeFin_nbr...
    elements_All.('Beam14_elements').('Element1').nodeIn_nbr,...
    elements_All.('Beam14_elements').(['Element' num2str(N_elem)]).nodeFin_nbr...
    ];
end

fprintf("\nModel reduction to %i DOFs \nCraig-Bampton's method",DOF_r_wanted);

% Extracting the DOFs corresponding to the nodes remaining
% Thanks to the exhaustive nodes list and the nodes numbers
DOF_r = nodes_All(nR_nodes_nbr(:),5:10);
% Sorting the array
DOF_r = sort(DOF_r(:));

% Initialising the condensed degrees of freedom
DOF_c = 1:nbr_DOFs;

% Suppressing from the full list the retained DOFs
for i = 1:length(DOF_r)
    for j = 1:length(DOF_c)
       if DOF_c(j) == DOF_r(i)
            DOF_c(j) = [];
            break;
       end
   end
end

Krr = K(DOF_r,DOF_r);
Krc = K(DOF_r,DOF_c);
Kcc = K(DOF_c,DOF_c);
Kcr = K(DOF_c,DOF_r);

Mrr = M(DOF_r,DOF_r);
Mcc = M(DOF_c,DOF_c);
Mrc = M(DOF_r,DOF_c);
Mcr = M(DOF_c,DOF_r);

Crr = C(DOF_r,DOF_r);
Ccc = C(DOF_c,DOF_c);
Crc = C(DOF_r,DOF_c);
Ccr = C(DOF_c,DOF_r);

Rcr = -Kcc\Kcr;

% Number of truncated modes to take into account for Craig Bampton method
tr_CB = 100;
fprintf(' - %i modes\n',tr_CB);
[Vcc,~] = eigs(Kcc,Mcc,tr_CB, 'sm');


R = [...
    eye(length(DOF_r)) , zeros(length(DOF_r),tr_CB); ...
    Rcr                ,                 Vcc         ...
    ];


K = [           ...
    Krr, Krc;   ...
    Kcr, Kcc;   ...
    ];

M = [           ...
    Mrr, Mrc;   ...
    Mcr, Mcc;   ...
    ];

C = [           ...
    Crr, Crc;   ...
    Ccr, Ccc;   ...
    ];



K_red = R' * K * R;
M_red = R' * M * R;
C_red = R' * C * R;


% Initialisation of the reduced external shaker force signal
p_red = p(DOF_r,:);
% Extension of the size to other modes taken into account
p_red(end+1:tr_CB + DOF_r_wanted,:) = 0;


t_Elapsed = toc(t_CB);
fprintf('%fs elapsed\n',t_Elapsed);

end