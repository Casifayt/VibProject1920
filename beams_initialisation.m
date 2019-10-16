%% Function initialising all beams elements
% It gives back a matrix (list of lists) beams which elements are the beams
% Each beams (row in the matrix) has seven elements :
%   - Initial node position
%   - Final node position
%   - Length [m]
%   - Orientation (unitary vector)
%   - Width [mm] (y-dimension)
%   - Thickness [mm] (x-dimension)
%   - Area[mm^2]

% The beams and nodes are ranked according to an increasing index of
% coordinates, taking x_1 = X, x_2 = Y, x_3 = Z

function [beams] = beams_initialisation(print)
% INPUTS
%   - print is a binary value to decide if the function prints its
%   iterations
% OUTPUTS
%   - beams is a list containing all beams with their properties

beams_made = 0;                         % Number of beams initialised
squareSize = 70;                        % size of the square section of beams

%% Beams in the (x,y,0) plane
% There are 13 beams in this plane

% Beams // 1_x
% There are 5 beams parallel to the x axis in the (x,y,0) plane
nbr_beams_this_case = 0;                

length = 4;                             % All those beams are 4m long
if print == 1
    fprintf('\nBEAMS // 1_x IN THE (x,y,0) PLANE\n');
end
for i = 0:4
   nodeIn = [0,5*i,0];
   nodeFin = [4,5*i,0];
   if print == 1
       fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   end
   orientation = [1,0,0];
   current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
   str = ['Beam' num2str(beams_made+(i+1))];
   beams.(str)=current_beam;
   %beams{beams_made+(i+1)}={nodeIn,nodeFin,length,orientation,};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;



% Beams // 1_y
% There are 8 beams parallel to the y axis in the (x,y,0) plane

length = 5;                             % All those beams are 5m long
if print == 1
    fprintf('\nBEAMS // 1_y IN (x,y,0) PLANE\n');
end
for i = 0:7
   if i < 4
       nodeIn = [0,5*i,0];
       nodeFin = [0,5*(i+1),0];
   else
       nodeIn = [4,5*(i-4),0];
       nodeFin = [4,5*(i-4+1),0];
   end
   if print == 1
       fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   end
   orientation = [0,1,0];
   current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
   str = ['Beam' num2str(beams_made+(i+1))];
   beams.(str)=current_beam;
   %beams{beams_made+(1+i)} = {nodeIn,nodeFin,length,orientation};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

%% Beams in the (x,y,3) plane
% There are 7 beams in this plane

% Beams // 1_x
% There are 3 beams parallel to the x axis in the (x,y,3) plane


length = 4;                             % All those beams are 4m long
if print == 1
    fprintf('\nBEAMS // 1_x IN THE (x,y,3) PLANE\n');
end
for i = 1:3
   nodeIn = [0,5*i,3];
   nodeFin = [4,5*i,3];
   if print == 1
       fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',i+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   end
   orientation = [1,0,0];
   current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
   str = ['Beam' num2str(beams_made+i)];
   beams.(str)=current_beam;
   %beams{i+beams_made}={nodeIn,nodeFin,length,orientation};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

% Beams // 1_y
% There are 4 beams parallel to the y axis in the (x,y,3) plane

length = 5;                             % All those beams are 5m long
if print == 1
    fprintf('\nBEAMS // 1_y IN THE (x,y,3) PLANE\n');
end
for i = 1:4
    if i < 3
       nodeIn = [0,5*i,3];
       nodeFin = [0,5*(i+1),3];
   else
       nodeIn = [4,5*((i-3)+1),3];
       nodeFin = [4,5*((i-3)+2),3];
    end
   if print == 1
       fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',i+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   end
   orientation = [0,1,0];
   current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
   str = ['Beam' num2str(beams_made+i)];
   beams.(str)=current_beam;
   %beams{i+beams_made} = {nodeIn,nodeFin,length,orientation};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

%% Beams in the (0,y,z) plane
% There are 7 beams in this plane

% Beams // 1_z axis
% There are 3 beams parallel to the z axis in the (0,y,z) plane

length = 3;                             % All those beams are 3m long
if print == 1
    fprintf('\nBEAMS // 1_z IN THE (0,y,z) PLANE\n');
end
for i = 1:3
    nodeIn = [0,5*i,0];
    nodeFin = [0,5*i,3];
    if print == 1
        fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',i+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    end
    orientation = [0,0,1];
    current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
    str = ['Beam' num2str(beams_made+i)];
    beams.(str)=current_beam;
    %beams{i+beams_made}={nodeIn,nodeFin,length,orientation};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

% Diagonal beams
% There are 4 diagonal beams in the (0,y,z) plane

length = sqrt(5^2+3^2);
if print == 1
    fprintf('\nBEAMS DIAGONAL IN THE (0,y,z) PLANE\n');
end
for i = 0:3
    if i == 0
       nodeIn = [0,0,0];
       nodeFin = [0,5,3];
       orientation = [0,1,1];
    elseif mod(i,2)~=0
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1),nodeIn(2)+5,nodeIn(3)-3];
       orientation = [0,1,-1];
    else
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1),nodeIn(2)+5,nodeIn(3)+3];
       orientation = [0,1,1];
    end
    if print == 1
        fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    end
    current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
    str = ['Beam' num2str(beams_made+(i+1))];
    beams.(str)=current_beam;
    %beams{(1+i)+beams_made}={nodeIn,nodeFin,length,orientation};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

%% Beams in the (4,y,z) plane
% There are 7 beams in this plane

% Beams // 1_z axis
% There are 3 beams parallel to the z axis in the (4,y,z) plane

length = 3;                             % All those beams are 3m long
if print == 1
    fprintf('\nBEAMS // 1_z IN THE (4,y,z) PLANE\n');
end
for i = 1:3
    nodeIn = [4,5*i,0];
    nodeFin = [4,5*i,3];
    if print == 1
        fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',i+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    end
    orientation = [0,0,1];
    current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
    str = ['Beam' num2str(beams_made+i)];
    beams.(str)=current_beam;
    %beams{i+beams_made}={nodeIn,nodeFin,length,orientation};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

% Diagonal beams
% There are 4 diagonal beams in the (4,y,z) plane

length = sqrt(5^2+3^2);
if print == 1
    fprintf('\nDIAGONAL BEAMS IN (4,y,z) plane\n');
end
for i = 0:3
    if i == 0
       nodeIn = [4,0,0];
       nodeFin = [4,5,3];
       orientation = [0,1,1];
    elseif mod(i,2)~=0
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1),nodeIn(2)+5,nodeIn(3)-3];
       orientation = [0,1,-1];
    else
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1),nodeIn(2)+5,nodeIn(3)+3];
       orientation = [0,1,1];
    end
    if print == 1
        fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    end
    current_beam = struct('node_Initial', nodeIn, 'node_Final', nodeFin, 'Length', length, 'Orientation', orientation);
    str = ['Beam' num2str(beams_made+(i+1))];
    beams.(str)=current_beam;
    %beams{(1+i)+beams_made}={nodeIn,nodeFin,length,orientation};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;
end