%% Function of discretisation of the beam in N elements
% This function takes a beam and returns a structure with 
% the N discretised elements of the given beam

function [elements] = discretisation(beam,N,beam_nbr,print)
% INPUTS
%  - beam       : Beam being discretised                    (structure)[/]
%  - N          : Total number of elements wanted           (int)[/]
%  - beam_nbr   : Index of the beam (for the prints)        (int)[/]
%  - print      : Binary value to print details or not      (int)[/]


% OUTPUTS
%  - elements   : All the discretised elements of the beam  (structure)[/]

%% Extraction of the beam properties
nodeInBeam = beam.node_Initial;     % Initial node
nodeFinBeam = beam.node_Final;      % Final node
lengthBeam = beam.Length;           % Length
orientation = beam.Orientation;     % Orientation
height = beam.Height;               % Height
width = beam.Width;                 % Width 


lengthEl = lengthBeam/N;            % Elements length       (float)[m]

if print == 1
    fprintf('\nDiscretisation of beam number %i into %i elements\n',beam_nbr,N);
end

for i = 1:N
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
    if print == 1
        fprintf('Element %i of length %f : nodeIn = (%i,%i,%i) and nodeFin = (%i,%i,%i)\n',...
            i,lengthEl, nodeInEl(1),nodeInEl(2),nodeInEl(3),nodeFinEl(1),nodeFinEl(2),nodeFinEl(3));
    end
    current_element = struct('node_Initial',nodeInEl, 'node_Final', nodeFinEl,'Length',...
        lengthEl,'Orientation',orientation,'Height',height,'Width',width);
    elements.(['Element' num2str(i)]) = current_element;
end
end