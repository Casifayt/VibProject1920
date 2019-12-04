%% Function constructing the cell array containing all nodes properties
% This function uses the properties stored in the fields of the structure
% elements_All to construct the exhaustive array of all nodes coordinates,
% and degrees of freedom

% Each node is stored in a 10 items list :
% x, y, z, 0, DOFx, DOFy, DOFz, DOFpsix, DOFpsiy, DOFpsiz
% Where (x,y,z) is the 3D location of the node in structural axis
% (DOFx, DOFy, DOFz) are the 3 translational DOFs in structural axis
% (DOFpsix, DOFpsiy, DOFpsiz) are the 3 rotational DOFs structural axis

% INPUTS
%   - elements_All  : structure with all the elements       (structure)[/]
%   - N_dof : number of total DOFs of the system            (int)[/]

% OUTPUT
%   - nodes_All : cell array with all the nodes properties  (cell array)[/]

function [nodes_All] = nodes_list_construction(elements_All,N_dof)

unit_vec_x = [1;0;0];       % Unit-vector needed to select elements parallel to x   (list)[/]
nodes_made = 1;             % Index of the number of nodes made so far              (int)[/]
dofs_made=1;

% Extraction of number of beams
num_beams = numel(fieldnames(elements_All)); 

% Extraction of number of elements on each beam
num_el = numel(fieldnames(elements_All.Beam1_elements));

% Initialisation of the nodes_All list
nodes_All = zeros(34*num_el-18,3 + 1 + N_dof);


for i = 1: num_beams                % Scanning on each beam
    current_beam = elements_All.(['Beam' num2str(i) '_elements']);  % Selecting one beam
    if num_el > 1                   % If we have discretised beams
        for j = 1:num_el       % Scanning on each element of this beam
            current_element = current_beam.(['Element' num2str(j)]);    % Selecting one element
            if current_element.Orientation*unit_vec_x == 1
            % If this element is parallel to x
                
                if j == 1
                %% If we are currently analysing first element of the beam
            
                %% Initial node
                nodes_All(nodes_made,1:3) = current_element.node_Initial(:);
                % Putting the coordinates of the initial node in the list
                
                % Constructing the DOFs of initial node :
                   for k = 1:N_dof
                      % Checking the constraint type of the node
                      if strcmp(current_element.nodeIn_cdt,'supported') == 1
                      % If supported then the z-dof is null 
                      % (avoided by k ~= 3)
                         if k ~= 3   
                             nodes_All(nodes_made,4+k) = dofs_made;
                             dofs_made=dofs_made+1;
                         end
                      elseif strcmp(current_element.nodeIn_cdt,'clamped') == 1
                      % If clamped then all translational dofs are null
                      % (avoided by k > 3)
                          if k > 3
                              nodes_All(nodes_made,4+k) = dofs_made;
                              dofs_made=dofs_made+1;
                          end
                      else
                      % If free then all dofs are directly good to go
                          nodes_All(nodes_made,4+k) = dofs_made;
                          dofs_made=dofs_made+1;
                      end
                   end
                   nodes_made = nodes_made+1;
                   
                   %% Final node
                   
                   nodes_All(nodes_made,1:3) = current_element.node_Final(:);
                   % Putting the coordinates of the final node in the list
                   
                   % Constructing the DOFs of initial node :
                   for k = 1:N_dof
                       % Checking the constraint type of the node
                       if strcmp(current_element.nodeFin_cdt,'supported') == 1
                       % If supported then the z-dof is null 
                       % (avoided by k ~= 3)    
                           if k ~= 3
                               nodes_All(nodes_made,4+k) = dofs_made;
                               dofs_made=dofs_made+1;
                           end
                       elseif strcmp(current_element.nodeFin_cdt,'clamped') == 1
                       % If clamped then all translational dofs are null
                       % (avoided by k > 3)    
                           if k > 3
                               nodes_All(nodes_made,4+k) = dofs_made;
                               dofs_made=dofs_made+1;
                           end
                       else
                       % If free then all dofs are directly good to go
                           nodes_All(nodes_made,4+k) = dofs_made;
                           dofs_made=dofs_made+1; 
                       end
                   end
                   nodes_made=nodes_made+1;
                else
                %% If we are not analysing the first element of the beam
                   nodes_All(nodes_made,1:3) = current_element.node_Final(:);     
                   for k = 1:N_dof
                       if strcmp(current_element.nodeFin_cdt,'supported') == 1
                           if k ~= 3
                               nodes_All(nodes_made,4+k) = dofs_made;
                               dofs_made=dofs_made+1;
                           end
                       elseif strcmp(current_element.nodeFin_cdt,'clamped') == 1
                           if k > 3
                               nodes_All(nodes_made,4+k) = dofs_made;
                               dofs_made=dofs_made+1;
                           end
                       else
                           nodes_All(nodes_made,4+k) = dofs_made;
                           dofs_made=dofs_made+1; 
                       end
                    end
                    nodes_made = nodes_made+1;
                end
            else
                if j ~= num_el
                   nodes_All(nodes_made,1:3) = current_element.node_Final(:);
                   for k = 1:N_dof
                       if strcmp(current_element.nodeFin_cdt,'supported') == 1
                           if k ~= 3
                               nodes_All(nodes_made,4+k) = dofs_made;
                               dofs_made=dofs_made+1;
                           end
                       elseif strcmp(current_element.nodeFin_cdt,'clamped') == 1
                           if k > 3
                               nodes_All(nodes_made,4+k) = dofs_made;
                               dofs_made=dofs_made+1;
                           end
                       else
                           nodes_All(nodes_made,4+k) = dofs_made;
                           dofs_made=dofs_made+1; 
                       end
                   end
                   nodes_made = nodes_made+1;
                end
            end
        end
    else                        % If we have only one element (the beam)
        current_element = current_beam.(['Element' num2str(num_el)]);    % Selecting one element
        if current_element.Orientation*unit_vec_x == 1
            
            nodes_All(nodes_made,1:3) = current_element.node_Initial(:);
            for k = 1:N_dof
               if strcmp(current_element.nodeIn_cdt,'supported') == 1
                   if k ~= 3
                       nodes_All(nodes_made,4+k) = dofs_made;
                       dofs_made=dofs_made+1;
                   end
               elseif strcmp(current_element.nodeIn_cdt,'clamped') == 1
                   if k > 3
                       nodes_All(nodes_made,4+k) = dofs_made;
                       dofs_made=dofs_made+1;
                   end
               else
                   nodes_All(nodes_made,4+k) = dofs_made;
                   dofs_made=dofs_made+1; 
               end
            end
            nodes_made = nodes_made+1;
            
            nodes_All(nodes_made,1:3) = current_element.node_Final(:);
            for k = 1:N_dof
               if strcmp(current_element.nodeIn_cdt,'supported') == 1
                   if k ~= 3
                       nodes_All(nodes_made,4+k) = dofs_made;
                       dofs_made=dofs_made+1;
                   end
               elseif strcmp(current_element.nodeIn_cdt,'clamped') == 1
                   if k > 3
                       nodes_All(nodes_made,4+k) = dofs_made;
                       dofs_made=dofs_made+1;
                   end
               else
                   nodes_All(nodes_made,4+k) = dofs_made;
                   dofs_made=dofs_made+1; 
               end
            end
            nodes_made=nodes_made+1;     
        end
    end
end
end