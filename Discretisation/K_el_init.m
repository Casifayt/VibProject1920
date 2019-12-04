%% Function initialising the elementary stiffness matrix
% This function initialises the stiffness mass matrix by taking
% advantage of the fields of the structure element given

function [K] = K_el_init(element,mat_prop)
% INPUTS
%   - element   : element to extract the matrix from        (structure)
%   - E         : Young's modulus                           (int)[Pa]
%   - G         : Shear modulus                             (float)[Pa]

% OUTPUT
%   - K         : Elementary stiffness matrix of the element (cell array)

%% Extraction of the element geometric properties
l = element.Length;         % Length                                    (int)[m]
h = element.Height*1e-3;    % Height                                    (int)[m]
b = element.Width*1e-3;     % Width                                     (int)[m]
A = h*b;                    % Computing cross section                   (float)[m^2]

%% Extraction of the material properties
E = mat_prop.E;             % Young's modulus                           (int)[Pa]
G = mat_prop.G;             % Shear modulus                             (float)[Pa]


if h==b     
    % SQUARE BEAMS
    J = 9*h^4/64;           % Torsional constant for square section     (float)[m^4]
    Iy = h^4/12;            % Second moment of area about local axis y  (float)[m^4]
    Iz = Iy;                % Second moment of area about local axis z  (float)[m^4]
else        
    % RECTANGULAR BEAM
 
    % Torsional constant for rectangular section                        (float)[m^4]
    % https://en.wikipedia.org/wiki/Torsion_constant
    J = h*b^3*(1/3 - 0.21*(b/h)*(1-(b^4/12/h^4))); 
    
    Iy = h*b^3/12;          % Second moment of area about local axis y  (float)[m^4]
    Iz = h^3*b/12;          % Second moment of area about local axis z  (float)[m^4]
end

%% Initialisation of the elementary stiffness matrix

y = E*Iy;             % Product of young's modulus and second
z = E*Iz;             % moment of area to save space

K = ...
[
E*A/l,          0,           0,      0,        0,         0, -E*A/l,          0,         0,      0,        0,        0; 
     0,  12*z/l^3,           0,      0,        0,   6*z/l^2,      0,  -12*z/l^3,         0,      0,        0,  6*z/l^2;
     0,         0,    12*y/l^3,      0, -6*y/l^2,         0,      0,          0, -12*y/l^3,      0, -6*y/l^2,        0;
     0,         0,           0,  G*J/l,        0,         0,      0,          0,         0, -G*J/l,        0,        0; 
     0,         0,    -6*y/l^2,      0,    4*y/l,         0,      0,          0,   6*y/l^2,      0,    2*y/l,        0;
     0,   6*z/l^2,           0,      0,        0,     4*z/l,      0,   -6*z/l^2,         0,      0,        0,    2*z/l;
-E*A/l,         0,           0,      0,        0,         0,  E*A/l,          0,         0,      0,        0,        0;
     0, -12*z/l^3,           0,      0,        0,  -6*z/l^2,      0,   12*z/l^3,         0,      0,        0, -6*z/l^2;
     0,         0,   -12*y/l^3,      0,  6*y/l^2,         0,      0,          0,  12*y/l^3,      0,  6*y/l^2,        0; 
     0,         0,           0, -G*J/l,        0,         0,      0,          0,         0,  G*J/l,        0,        0;
     0,         0,    -6*y/l^2,      0,    2*y/l,         0,      0,          0,   6*y/l^2,      0,    4*y/l,        0;
     0,   6*z/l^2,           0,      0,        0,     2*z/l,      0,   -6*z/l^2,         0,      0,        0,    4*z/l;
];
end