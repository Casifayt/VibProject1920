%% Function initialising the elementary stiffness matrix
% This function initialises the stiffness mass matrix by taking
% advantage of the fields of the structure element given

function [K] = K_el_init(element,E,G)
% INPUTS
%   - element   : element to extract the matrix from        (structure)
%   - E         : Young's modulus                           (int)[Pa]
%   - G         : Shear modulus                             (float)[Pa]

% OUTPUT
%   - K         : Elementary stiffness matrix of the element (cell array)

%% Extraction of the element properties
l = element.Length;         % Length (l)
h = element.Height;         % Height (mm)
b = element.Width;          % Width (mm)

%% Computation of the element properties needed
A = h*b;                    % Computing cross section are (mm^2)

if h==b     
    % SQUARE BEAMS
    J = 9*h^4/64;           % Torsional constant for square section (mm^4)
    Iy = h^4/12;            % Second moment of area about local axis y (mm^4)
    Iz = Iy;                % Second moment of area about local axis z (mm^4)
else        
    % RECTANGULAR BEAM
 
    % Torsional constant for rectangular section (mm^4)
    % https://en.wikipedia.org/wiki/Torsion_constant
    J = h*b^3*(1/3 - 0.21*(b/h)*(1-(b^4/12/h^4))); 
    
    Iy = h*b^3/12;          % Second moment of area about local axis y (mm^4)
    Iz = h^3*b/12;          % Second moment of area about local axis z (mm^4)
end

%% Initialisation of the elementary stiffness matrix

y = E*Iy*1e-12;             % Rewriting of young's modulus and second moment of area
z = E*Iz*1e-12;             % to save space and computation factor to N.m^2

K = ...
[
E*A/l,          0,           0,      0,        0,         0, -E*A/l,          0,         0,      0,        0,        0; 
     0,  12*z/l^3,           0,      0,        0,  6*z/l^2,       0, -12*z/l^3,          0,      0,        0,  6*z/l^2;
     0,         0,    12*y/l^3,      0, -6*y/l^2,         0,      0,          0, -12*y/l^3,      0, -6*y/l^2,        0;
     0,         0,           0,  G*J/l,        0,         0,      0,          0,         0, -G*J/l,        0,        0; 
     0,         0,    -6*y/l^2,      0,    4*y/l,         0,      0,          0,   6*y/l^2,      0,    2*y/l,        0;
     0,   6*z/l^2,           0,      0,        0,    4*z/l,       0,  -6*z/l^2,          0,      0,        0,    2*z/l;
-E*A/l,         0,           0,      0,        0,         0,  E*A/l,          0,         0,      0,        0,        0;
     0, -12*z/l^3,           0,      0,        0, -6*z/l^2,       0,  12*z/l^3,          0,      0,        0, -6*z/l^2;
     0,         0,   -12*y/l^3,      0,  6*y/l^2,         0,      0,          0,  12*y/l^3,      0,  6*y/l^2,        0; 
     0,         0,           0, -G*J/l,        0,         0,      0,          0,         0,  G*J/l,        0,        0;
     0,         0,    -6*y/l^2,      0,    2*y/l,         0,      0,          0,   6*y/l^2,      0,    4*y/l,        0;
     0,   6*z/l^2,           0,      0,        0,    2*z/l,       0,  -6*z/l^2,          0,      0,        0,    4*z/l;
];
end