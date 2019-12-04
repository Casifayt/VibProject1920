%% Function constructing the Rotation matrices
% This function computes all the 3D 6-DOFs rotation matrices needed 
% in the current structure. It produces a structure which fields are :
%   - Rx : rotation matrix for x-oriented elements
%   - Ry : rotation matrix for y-oriented elements
%   - Rz : rotation matrix for z-oriented elements
%   - Rdiag_up : rotation matrix for diagonal up-oriented elements
%   - Rdiag_down : rotation matrix for diagonal down-oriented elements

% INPUTS : /

% OUTPUT
%   - rot_matrices :    structure with all the rotation matrices            (structure)[/]

function [rot_matrices] = rot_mat_init()

angle = atan(3/5);
CA = cos(angle);
SA = sin(angle);

% Cosinus and sinus of pi/2
CP = 0; 
SP = 1;

% Cosinus and sinus of pi
CPP = -1;
SPP = 0;


%% X axis
rot_mat_x = [
                1,   0,   0;
                0,  CP,  SP;
                0, -SP,  CP;
            ];

Rx = [
 rot_mat_x, zeros(3,9);
 zeros(3,3), rot_mat_x, zeros(3,6);
 zeros(3,6), rot_mat_x, zeros(3,3);
 zeros(3,9), rot_mat_x;
 ];


%% Y axis
Ry1 = Rx;

rot_mat_y2 = [
                CP,  SP,  0;
               -SP,  CP,  0;
                 0,  0,   1;
            ];
        
Ry2 = [
 rot_mat_y2, zeros(3,9);
 zeros(3,3), rot_mat_y2, zeros(3,6);
 zeros(3,6), rot_mat_y2, zeros(3,3);
 zeros(3,9), rot_mat_y2;
 ];
Ry = Ry1*Ry2;


%% Z axis
rot_mat_z1 = [
                1,    0,    0;
                0,  CPP,  SPP;
                0, -SPP,  CPP;
            ];

Rz1 = [
 rot_mat_z1, zeros(3,9);
 zeros(3,3), rot_mat_z1, zeros(3,6);
 zeros(3,6), rot_mat_z1, zeros(3,3);
 zeros(3,9), rot_mat_z1;
 ];

rot_mat_z2 = [
                CP,  0, -SP;
                 0,  1,   0;
                SP,  0,  CP;
            ];
        
Rz2 = [
 rot_mat_z2, zeros(3,9);
 zeros(3,3), rot_mat_z2, zeros(3,6);
 zeros(3,6), rot_mat_z2, zeros(3,3);
 zeros(3,9), rot_mat_z2;
 ];

Rz = Rz1 * Rz2;


%% Diag up
rot_mat_diag_up_1 = [ 
                    CA, SA, 0;
                   -SA, CA, 0;
                     0, 0,  1;
                    ];

                
R_diag_up1 = [
 rot_mat_diag_up_1, zeros(3,9);
 zeros(3,3), rot_mat_diag_up_1, zeros(3,6);
 zeros(3,6), rot_mat_diag_up_1, zeros(3,3);
 zeros(3,9), rot_mat_diag_up_1;
 ];


R_diag2 = Rx;
R_diag3 = Ry2;

R_diag_up = R_diag_up1 * R_diag2 * R_diag3;

%% Diag down

rot_mat_diag_up_1 = [ 
                    CA,-SA, 0;
                    SA, CA, 0;
                     0, 0,  1;
                    ];
R_diag_down1 = [
 rot_mat_diag_up_1, zeros(3,9);
 zeros(3,3), rot_mat_diag_up_1, zeros(3,6);
 zeros(3,6), rot_mat_diag_up_1, zeros(3,3);
 zeros(3,9), rot_mat_diag_up_1;
 ];

R_diag_down = R_diag_down1 * R_diag2 * R_diag3;

rot_matrices =  struct('Rx', Rx, 'Ry', Ry, 'Rz', Rz, 'Rdiag_up', R_diag_up, 'Rdiag_down', R_diag_down);
end