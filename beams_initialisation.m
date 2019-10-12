%% Function initialising all beams elements
% It gives back a matrix (list of lists) beams which elements are the beams
% Each beams (row in the matrix) has three elements :
%   - Initial node
%   - Final node
%   - Length
% The beams and nodes are ranked according to an increasing index of
% coordinates, taking x_1 = X, x_2 = Y, x_3 = Z

function [beams] = beams_initialisation()

beams_made = 0;                         % Number of beams initialised

%% Beams in the (x,y,0) plane
% There are 13 beams in this plane

% Beams // 1_x
% There are 5 beams parallel to the x axis in the (x,y,0) plane
nbr_beams_this_case = 0;                

length = 4;                             % All those beams are 4m long
fprintf('\nBEAMS // 1_x IN THE (x,y,0) PLANE\n');
for i = 0:4
   nodeIn = [0;5*i;0];
   nodeFin = [4;5*i;0];
   fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   beams{beams_made+(i+1)}={nodeIn,nodeFin,length};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;



% Beams // 1_y
% There are 8 beams parallel to the y axis in the (x,y,0) plane

length = 5;                             % All those beams are 5m long
fprintf('\nBEAMS // 1_y IN (x,y,0) PLANE\n');
for i = 0:7
   if i < 4
       nodeIn = [0;5*i;0];
       nodeFin = [0;5*(i+1);0];
   else
       nodeIn = [4;5*(i-4);0];
       nodeFin = [4;5*(i-4+1);0];
   end
   fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   beams{beams_made+(1+i)} = {nodeIn,nodeFin,length};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

%% Beams in the (x,y,3) plane
% There are 7 beams in this plane

% Beams // 1_x
% There are 3 beams parallel to the x axis in the (x,y,3) plane


length = 4;                             % All those beams are 4m long
fprintf('\nBEAMS // 1_x IN THE (x,y,3) PLANE\n');
for i = 1:3
   nodeIn = [0;5*i;3];
   nodeFin = [4;5*i;3];
   fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
   beams{i+beams_made}={nodeIn,nodeFin,length};
   nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

% Beams // 1_y
% There are 4 beams parallel to the y axis in the (x,y,3) plane

length = 5;                             % All those beams are 5m long
fprintf('\nBEAMS // 1_y IN THE (x,y,3) PLANE\n');
for i = 1:4
    if i < 3
       nodeIn = [0;5*i;3];
       nodeFin = [0;5*(i+1);3];
   else
       nodeIn = [4;5*((i-3)+1);3];
       nodeFin = [4;5*((i-3)+2);3];
    end
   fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    beams{i+beams_made} = {nodeIn,nodeFin,length};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

%% Beams in the (0,y,z) plane
% There are 7 beams in this plane

% Beams // 1_z axis
% There are 3 beams parallel to the z axis in the (0,y,z) plane

length = 3;                             % All those beams are 3m long
fprintf('\nBEAMS // 1_z IN THE (0,y,z) PLANE\n');
for i = 1:3
    nodeIn = [0;5*i;0];
    nodeFin = [0;5*i;3];
    fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    beams{i+beams_made}={nodeIn,nodeFin,length};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

% Diagonal beams
% There are 4 diagonal beams in the (0,y,z) plane

length = sqrt(5^2+3^2);
fprintf('\nBEAMS DIAGONAL IN THE (0,y,z) PLANE\n');
for i = 0:3
    if i == 0
       nodeIn = [0;0;0];
       nodeFin = [0;5;3];
    elseif mod(i,2)~=0
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1);nodeIn(2)+5;nodeIn(3)-3];
    else
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1);nodeIn(2)+5;nodeIn(3)+3];
    end
    fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    beams{(1+i)+beams_made}={nodeIn,nodeFin,length};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

%% Beams in the (4,y,z) plane
% There are 7 beams in this plane

% Beams // 1_z axis
% There are 3 beams parallel to the z axis in the (4,y,z) plane

length = 3;                             % All those beams are 3m long
fprintf('\nBEAMS // 1_z IN THE (4,y,z) PLANE\n');
for i = 1:3
    nodeIn = [4;5*i;0];
    nodeFin = [4;5*i;3];
    fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    beams{i+beams_made}={nodeIn,nodeFin,length};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

% Diagonal beams
% There are 4 diagonal beams in the (4,y,z) plane

length = sqrt(5^2+3^2);
fprintf('\nDIAGONAL BEAMS IN (4,y,z) plane\n');
for i = 0:3
    if i == 0
       nodeIn = [4;0;0];
       nodeFin = [4;5;3];
    elseif mod(i,2)~=0
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1);nodeIn(2)+5;nodeIn(3)-3];
    else
       nodeIn = nodeFin;
       nodeFin = [nodeIn(1);nodeIn(2)+5;nodeIn(3)+3];
    end
    fprintf('Beams %i with initial node (%i,%i,%i) and final node (%i,%i,%i)\n',(1+i)+beams_made,nodeIn(1),nodeIn(2),nodeIn(3),nodeFin(1),nodeFin(2),nodeFin(3));
    beams{(1+i)+beams_made}={nodeIn,nodeFin,length};
    nbr_beams_this_case = nbr_beams_this_case + 1;
end
beams_made = beams_made + nbr_beams_this_case;
nbr_beams_this_case = 0;

end