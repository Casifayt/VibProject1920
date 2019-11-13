%% Function of discretisation of the beam in N elements
% This function takes the structure with all the beams and returns
% a structure with the N discretised elements of the given beam. 

% Each element has the following fields : 
%   - Initial node position         (list)[/]
%   - Initial node condition        (string)[/]
%   - Final node position           (list)[/]
%   - Final node condition          (string)[/]
%   - Length                        (int)[m]
%   - Orientation                   (list)[/]
%   - Height, local z-dimension     (int)[mm]
%   - Width, local x-dimension      (int)[mm]

% INPUTS
%  - beam_All   : Structure containing all the beams        (structure)[/]
%  - N_elem     : Number of elements wanted per beam        (int)[/]
%  - print      : Binary value to print details or not      (int)[/]


% OUTPUTS
%  - elements   : All the discretised elements of the beam  (structure)[/]

function [elements_All] = discretisation(beams_All,N_elem,print)

% Initialisation of the timer
t_discr = tic;

fprintf('\nDiscretisation into %i elements\n',N_elem);

nodes_made = 0;


for k = 1:numel(fieldnames(beams_All))
    beam = beams_All.(['Beam' num2str(k)]);

    %% Extraction of the beam properties
    nodeInBeam = beam.node_Initial;     % Initial node              (list)[/]
    nodeFinBeam = beam.node_Final;      % Final node                (list)[/]
    lengthBeam = beam.Length;           % Length                    (int)[m]
    orientation = beam.Orientation;     % Orientation               (list)[/]
    height = beam.Height;               % Height                    (int)[mm]
    width = beam.Width;                 % Width                     (int)[mm]
    
    
    lengthEl = lengthBeam/N_elem;            % Elements length           (float)[m]
    
    if print == 1
        fprintf('\nDiscretisation of beam number %i into %i elements\n',k,N_elem);
    end

        
    %% Discretisation
    for i = 1:N_elem
        nodeInEl = nodeInBeam; 
        nodeFinEl = nodeInBeam;
        if orientation == [1,0,0]       % Beams // 1_x
           nodeInEl(1) = nodeInBeam(1)+(i-1)*lengthEl;
           nodeFinEl(1) = nodeFinEl(1)+i*lengthEl; 
        elseif orientation == [0,1,0]    % Beams // 1_y
           nodeInEl(2) = nodeInBeam(2)+(i-1)*lengthEl;
           nodeFinEl(2) = nodeFinEl(2)+i*lengthEl;
        elseif orientation == [0,0,1]   % Beams // 1_z
           nodeInEl(3) = nodeInBeam(3)+(i-1)*lengthEl;
           nodeFinEl(3) = nodeFinEl(3)+i*lengthEl;
        elseif orientation == [0,1,1]   % Beams diagonal up
           nodeInEl(2) = nodeInBeam(2)+(i-1)*lengthEl*abs(nodeFinBeam(2)-nodeInBeam(2))/lengthBeam;
           nodeInEl(3) = nodeInBeam(3)+(i-1)*lengthEl*abs(nodeFinBeam(3)-nodeInBeam(3))/lengthBeam;
           nodeFinEl(2) = nodeInBeam(2)+i*lengthEl*abs(nodeFinBeam(2)-nodeInBeam(2))/lengthBeam;
           nodeFinEl(3) = nodeInBeam(3)+i*lengthEl*abs(nodeFinBeam(3)-nodeInBeam(3))/lengthBeam;
        elseif orientation == [0,1,-1]  % Beams diagonal down
           nodeInEl(2) = nodeInBeam(2)+(i-1)*lengthEl*abs(nodeFinBeam(2)-nodeInBeam(2))/lengthBeam;
           nodeInEl(3) = nodeInBeam(3)-(i-1)*lengthEl*abs(nodeFinBeam(3)-nodeInBeam(3))/lengthBeam;
           nodeFinEl(2) = nodeInBeam(2)+i*lengthEl*abs(nodeFinBeam(2)-nodeInBeam(2))/lengthBeam;
           nodeFinEl(3) = nodeInBeam(3)-i*lengthEl*abs(nodeFinBeam(3)-nodeInBeam(3))/lengthBeam;
        end

        % To avoid residus of order of the computer precision
        tol = 1e-12;                    % Tolerance on computations (float)[m]
        
        if nodeFinEl(2) < tol
           nodeFinEl(2) = 0;
        elseif nodeFinEl(3) < tol
           nodeFinEl(3) = 0;
        end

        % Implementation of nodes freedom conditions
        if N_elem == 1
            nodeIn_cdt = beam.nodeIn_cdt;
            nodeFin_cdt = beam.nodeFin_cdt;
        else    
            if i == 1
                nodeIn_cdt = beam.nodeIn_cdt;
                nodeFin_cdt = 'free';
            elseif i == N_elem
                nodeIn_cdt = 'free';
                nodeFin_cdt = beam.nodeFin_cdt;
            else
               nodeIn_cdt = 'free';
               nodeFin_cdt = 'free';
            end
        end
        
        if print == 1
            fprintf(['Element %i of length %f with initial node (%i,%i,%i) ' nodeIn_cdt ' and final node (%i,%i,%i) '...
                nodeFin_cdt '\n'],i,lengthEl,nodeInEl(1),nodeInEl(2),nodeInEl(3),nodeFinEl(1),nodeFinEl(2),nodeFinEl(3));
        end
        
        node_in_list = 0;
        
        if nodes_made ~=0
            for h = 1:nodes_made
                if abs(nodes_done_list(h,1) - nodeInEl(1)) < tol && abs(nodes_done_list(h,2) - nodeInEl(2)) < tol && abs(nodes_done_list(h,3) - nodeInEl(3)) < tol
                    node_in_list = h;
                    break;
                end
            end
            
            if node_in_list ~= 0
                nodeIn_nbr = node_in_list;
            else
                nodes_made = nodes_made + 1;
                nodes_done_list(nodes_made,:) = nodeInEl(:);
                nodeIn_nbr = nodes_made;
            end
        else
            nodes_done_list(1,:) = nodeInEl(:);
            nodes_made = nodes_made + 1;
            nodeIn_nbr = 1;
        end
        
        node_in_list = 0;
        
        for h = 1:nodes_made
            if abs(nodes_done_list(h,1) - nodeFinEl(1)) < tol && abs(nodes_done_list(h,2) - nodeFinEl(2)) < tol && abs(nodes_done_list(h,3) - nodeFinEl(3)) < tol
                node_in_list = h;
                break;
            end
        end

        if node_in_list ~= 0
            nodeFin_nbr = node_in_list;
        else
            nodes_made = nodes_made + 1;
            nodes_done_list(nodes_made,:) = nodeFinEl(:);
            nodeFin_nbr = nodes_made;
        end
        
        current_element = struct('node_Initial', nodeInEl, 'nodeIn_cdt',nodeIn_cdt, 'nodeIn_nbr', ...
            nodeIn_nbr, 'node_Final', nodeFinEl,'nodeFin_cdt', nodeFin_cdt,'nodeFin_nbr', nodeFin_nbr,...
            'Length',lengthEl, 'Orientation', orientation,'Height',height,'Width',width);
        elements.(['Element' num2str(i)]) = current_element;
    end
    elements_All.(['Beam' num2str(k) '_elements'])=elements;
end

tElapsed = toc(t_discr);
fprintf('(performed in %.2fs)\n',tElapsed);

end