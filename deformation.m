%% Function constructing the deformed elements according to the mode
% This function works as discretisation.m but uses the values calculated in
% the eigenvalue problem to compute the 3D-displacements of initial and final
% nodes of each element

% INPUTS
%   - mode : number of the desired mode                         (int)[/]
%   - V : structural displacements matrix                       (cell array)[/]
%   - locel : localisation element matrix                       (cell array)[/]
%   - elements_All : structure with all initial elements        (structure)[/]

% OUTPUT
%   - elements_Deformed : structure with all deformed elements  (structure)[/]

function [elements_Deformed] = deformation(mode, wColInV, V, locel, elements_All)

if mode > numel(wColInV)
  error('Choose a mode in the chosen truncated ones (mode <= tr).');
end


for i = 1:numel(fieldnames(elements_All))                   % Scanning on each beam
   
   beam = elements_All.(['Beam' num2str(i) '_elements']);   % Selecting one beam
   N_elem = numel(fieldnames(beam));                        % Extracting number of elements per beam
   
   for j = 1:N_elem                                         % Scanning on each element of the beam
       element = beam.(['Element' num2str(j)]);             % Selecting one element
       
       % Extraction of unchanged properties
       nodeIn_cdt = element.nodeIn_cdt;
       nodeFin_cdt = element.nodeFin_cdt;
       nodeIn_nbr = element.nodeIn_nbr;
       nodeFin_nbr = element.nodeFin_nbr;
       orientation = element.Orientation;
       height = element.Height;
       width = element.Width;
       
       
       el_nbr = (i-1)*N_elem + j;                           % Nbr of the current element
       
       
       % Computation of the deformations for the initial node according to its condition
       if strcmp(element.nodeIn_cdt, 'supported') == 1
           def_In_trans  = [V(locel(el_nbr,1:2),wColInV(mode)); 0];
       elseif strcmp(element.nodeIn_cdt, 'clamped') == 1
           def_In_trans = zeros(1,3);
       else
           def_In_trans  = V(locel(el_nbr,1:3),wColInV(mode));
       end
       
       % Computation of the deformations for the final node according to its condition
       if strcmp(element.nodeFin_cdt, 'supported') == 1
           def_Fin_trans  = [V(locel(el_nbr,7:8),wColInV(mode)); 0];
       elseif strcmp(element.nodeFin_cdt, 'clamped') == 1
           def_Fin_trans = zeros(1,3);
       else
           def_Fin_trans  = V(locel(el_nbr,7:9),wColInV(mode));
       end
       
       def_In_trans = def_In_trans';
       def_Fin_trans = def_Fin_trans';
       
       % Computation of nodes coordinates according to undeformed node and displacement
       node_Initial(:) = element.node_Initial(:) + def_In_trans(:);
       node_Final(:) = element.node_Final(:) + def_Fin_trans(:);
       
       % Computation of the length of the deformed element
       length = sqrt(...
           (node_Initial(1)-node_Final(1))^2 ...
       + (node_Initial(2)-node_Final(2))^2 ...
       + (node_Initial(3)-node_Final(3))^2);
       
       % Storing all informations into a structure
       current_element = struct('node_Initial', node_Initial, 'nodeIn_cdt',nodeIn_cdt, 'nodeIn_nbr', ...
            nodeIn_nbr, 'node_Final', node_Final,'nodeFin_cdt', nodeFin_cdt,'nodeFin_nbr', nodeFin_nbr,...
            'Length',length, 'Orientation', orientation,'Height',height,'Width',width);
        
        elements.(['Element' num2str(j)]) = current_element;
   end
   elements_Deformed.(['Beam' num2str(i) '_elements'])=elements;
end
end