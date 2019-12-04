%% Function constructing the structural stiffness and mass matrices
% This function scans each element, builds the corresponding elementary
% stiffness and mass matrices, and then constructs the structural stiffness
% and mass matrices according to the nodes condition(s) of this element

% INPUTS
%   - elements_All  : Structure containing all elements         (struct)[/]
%   - rot_matrices  : Structure with all rotation matrices      (struct)[/]
%   - mat_prop      : Structure containing material properties  (struct)[/]
%   - nodes_All     : Array containing all nodes and properties (array)[/]
%   - locel         : Localisation element matrix               (array)[/]

% OUTPUTS
%   - K_S           : Stiffness structural matrix               (array)[/]
%   - M_S           : Mass structural matrix                    (array)[/]

function [K_S,M_S] = struct_mat_init(elements_All, rot_matrices, mat_prop, nodes_All, locel)

N_elem = numel(fieldnames(elements_All.('Beam1_elements')));


%% Initialisation of structural matrices
total_nbr_dofs = nodes_All(end,end);        % Extraction of the total number of degrees of freedom
K_S = zeros(total_nbr_dofs,total_nbr_dofs); % Initialisation of the structural stiffness matrix
M_S = zeros(total_nbr_dofs,total_nbr_dofs); % Initialisation of the structural mass matrix


for i = 1:numel(fieldnames(elements_All))                           % Scanning on each beam
    current_beam = elements_All.(['Beam' num2str(i) '_elements']);  % Selecting one beam
    for j = 1:numel(fieldnames(current_beam))                       % Scanning on each element of this beam
        current_element = current_beam.(['Element' num2str(j)]);    % Selecting one element
        
        %% Construction of elementary matrices
        
        % Construction of the corresponding elementary matrices
        K_el = K_el_init(current_element,mat_prop);
        M_el = M_el_init(current_element,mat_prop);
        
        % Definition of the rotation matrix needed according to element orientation
        if current_element.Orientation(1) == 0
            if current_element.Orientation(2) == 1
                if current_element.Orientation(3) == 1
                  R = rot_matrices.Rdiag_up;
                elseif current_element.Orientation(3) == -1
                  R = rot_matrices.Rdiag_down;
                else
                  R = rot_matrices.Ry;
                end
            else
               R = rot_matrices.Rz;
            end
        else
            R = rot_matrices.Rx;
        end
        
        % Rotation of elementary matrices
        K_el = R' * K_el * R;
        M_el = R' * M_el * R;
        
        %% Construction of structural matrices
        el_nbr = (i-1)*N_elem + j;                  % Nbr of the current element
        
        
        % Taking into account initial node condition to build the index
        % of the locel matrix to be used to compute structural matrix value
        if strcmp(current_element.nodeIn_cdt,'supported') == 1 
            indexIn =  [1,2,4,5,6];
        elseif strcmp(current_element.nodeIn_cdt,'clamped') == 1 
            indexIn = 4:6;
        else
            indexIn = 1:6;
        end
        
        % Taking into account final node condition to build the index
        % of the locel matrix to be used to compute structural matrix value
        if strcmp(current_element.nodeFin_cdt,'supported') == 1 
            indexFin =  [7,8,10,11,12];
        elseif strcmp(current_element.nodeFin_cdt,'clamped') == 1 
            indexFin = 10:12;
        else
            indexFin = 7:12;
        end
        
        % Concatenation of index for initial node and index for final node
        indexTot = [indexIn indexFin];
         
        % Construction of structural matrix
        K_S(locel(el_nbr,indexTot),locel(el_nbr,indexTot))= K_S(locel(el_nbr,indexTot),locel(el_nbr,indexTot))  + K_el(indexTot,indexTot);
        M_S(locel(el_nbr,indexTot),locel(el_nbr,indexTot))= M_S(locel(el_nbr,indexTot),locel(el_nbr,indexTot))  + M_el(indexTot,indexTot);
    end
end
end