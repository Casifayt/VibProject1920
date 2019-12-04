%% Function of model redution using the static Guyan Irons' method
% This function reduces the problem size into a certain number of retained
% nodes. This reductions costs in accuracy but yields a high decrease in
% computation time.

% INPUTS
%   - elements_All : Structure with all elements properties     (structure)[/]
%   - nodes_All :    Array with all nodes and their DOFs        (structure)[/]
%   - K_S :          Structural stiffness matrix                (matrix)[/]
%   - M_S :          Structural mass matrix                     (matrix)[/]
%   - C :            Damping matrix                             (matrix)[/]
%   - V :            Matrix of eigenmodes                       (matrix)[/]
%   - wColInV :      Array of corresponding columns of modes in V(array)[/]
%   - tr :           Number of truncated modes                  (int)[/]
%   - p :            External excitation force vector           (matrix)[/]
%   - DOF_R_wanted : Number of retained degrees of freedom      (int)[/]

% OUTPUTS
%   - K_red :       Reduced stiffness matrix                    (matrix)[/]
%   - M_red :       Reduced mass matrix                         (matrix)[/]
%   - C_red :       Reduced damping matrix                      (matrix)[/]
%   - p_red :       Reduced external shaker force               (matrix)[/]

function [K_red,M_red,C_red,p_red] = guyan_irons(elements_All, nodes_All, K_S, M_S, C, ...
    V, wColInV, tr, p, DOF_r_wanted)

% Initialising timer
t_GI = tic;

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

nbr_nodes = length(nR_nodes_nbr);
nbr_DOFs_reduced = 6 *nbr_nodes;

fprintf("\nModel reduction to %i DOFs - Guyan-Irons' method\n",...
    nbr_DOFs_reduced);

% Extracting the DOFs corresponding to the nodes remaining
% Thanks to the exhaustive nodes list and the nodes numbers
DOF_r = nodes_All(nR_nodes_nbr,5:10);
% Sorting the array and making it a (1,24) array
DOF_r = sort(DOF_r(:));
DOF_r = DOF_r';


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



Krr = K_S(DOF_r,DOF_r);
Krc = K_S(DOF_r,DOF_c);
Kcc = K_S(DOF_c,DOF_c);
Kcr = K_S(DOF_c,DOF_r);


Rcr = - Kcc\Kcr;


R = [                   ...
    eye(length(DOF_r)); ...
    Rcr;                ...
    ];


Mrr = M_S(DOF_r,DOF_r);
Mcc = M_S(DOF_c,DOF_c);
Mrc = M_S(DOF_r,DOF_c);
Mcr = M_S(DOF_c,DOF_r);

Crr = C(DOF_r,DOF_r);
Ccc = C(DOF_c,DOF_c);
Crc = C(DOF_r,DOF_c);
Ccr = C(DOF_c,DOF_r);


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


% Extracting sorted eigenmode shapes
xR = R\V;
xC = -(Kcc\Kcr) * xR;

% Verification of the assumption F_C \approx 0
min_FR = zeros(1,tr);
max_FC = zeros(1,tr);
problem = 0;

for mode = 1:tr
     xR_temp = xR(:,wColInV(mode));
     xC_temp = xC(:,wColInV(mode));
    
    min_temp_FR = min(min(abs(Krr * xR_temp + Krc * xC_temp)));
    if min_temp_FR > min_FR(mode)
        min_FR(mode) = min_temp_FR;
    end
    
    max_temp_FC = max(max(abs(Kcr * xR_temp + Kcc * xC_temp)));
    if max_temp_FC > max_FC(mode)
        max_FC(mode) = max_temp_FC;
    end
    
    if min_FR(mode)/max_FC(mode) < 100
        fprintf('\nGuyan Irons assumption calculated wrong, at mode nbr %i\n',mode);
        problem = 1;
    end
end
if problem ~= 1
    fprintf('Minimal ratio Fr/Fc = %e\n',min(min_FR)/max(max_FC));
    fprintf('Maximal value of Fc = %e\n',max(max_FC));
end

% Reduction of external shaker force signal
p_red = p(DOF_r,:);


t_Elapsed = toc(t_GI);
fprintf('%fs elapsed\n',t_Elapsed);

end