%% Function constructing the localization element matrix
% This matrix gives the DOF corresponding to an element by defining
% array (e,j) as the degree of freedom q_{e,j} of element e
% and is needed to construct the structural stiffness and mass matrices

% INPUTS
%   - elements_All  : structure of all the elements             (structure)[/]
%   - nodes_All     : cell array containing all node properties (cell array)[/]

% OUTPUT
%   - locel         : localization element matrix               (matrix)[/]

function [locel] = locel_matrix_init(elements_All, nodes_All)

N_beam = numel(fieldnames(elements_All));
N_elem = numel(fieldnames(elements_All.('Beam1_elements')));

locel = zeros(N_elem,12);

for i = 1:N_beam               
    % Scanning each beam
    current_beam = elements_All.(['Beam' num2str(i) '_elements']);  % Selecting one beam
    
    for j = 1:N_elem           
        % Scanning each element of this beam
        current_element = current_beam.(['Element' num2str(j)]);    % Selecting one element
        el_nbr = (i-1)*N_elem + j;
        
        nodeIn_nbr = current_element.nodeIn_nbr;
        % Extracting the index number of the initial node of the element
        nodeFin_nbr = current_element.nodeFin_nbr;
        % Extracting the index number of the final node of the element
        
        % As the list nodes_All contains all node properties for node
        % number x at row x, presented as 
        % (x,y,z, DOFx, DOFy, DOFz, DOF_Psix, DOF_Psiy, DOF_Psi_z)
        % the DOFs of the node are given by
        
        locel(el_nbr,1:6)  = nodes_All(nodeIn_nbr,5:10);
        locel(el_nbr,7:12) = nodes_All(nodeFin_nbr,5:10);
    end
end
end