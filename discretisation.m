%% Function of discretisation of the beam in N elements
% This function takes a beam in the form of a cell array 
% as initialised by the beams_initialisation function
% It gives back a cell array of the N elements encrypted with the same
% routine
%   - Initial node
%   - Final node
%   - Length

function [elements] = discretisation(beam,N,beam_nbr,print)
% INPUTS
%  - beam is the beam element to be discretised in the form of a structure
%  - N is the number of elements in which the beam will be discretised
%  - beam_nbr is the number of the beam being discretised
%  - print is a str value to choose if the function will print its
%  iterations
%
% OUTPUTS
%  - elements is a cell array containing details about the elements of the
%  beam being discretised

nodeInBeam = beam.node_Initial;        % Extracting initial node from beam
nodeFinBeam = beam.node_Final;       % Extracting final node from beam
lengthBeam = beam.Length;        % Extracting beam length
orientation = beam.Orientation;       % Extracting beam orientation

lengthEl = lengthBeam/N;
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
        fprintf('Element %i of length %f : nodeIn = (%i,%i,%i) and nodeFin = (%i,%i,%i)\n',i,lengthEl, nodeInEl(1),nodeInEl(2),nodeInEl(3),nodeFinEl(1),nodeFinEl(2),nodeFinEl(3));
    end
    current_element = struct('node_Initial',nodeInEl, 'node_Final', nodeFinEl,'Length',lengthEl,'Orientation',orientation);
    elements.(['Element' num2str(i)]) = current_element;
end
end