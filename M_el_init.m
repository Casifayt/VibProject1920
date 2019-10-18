%% Function initialising the elementary mass matrix
% This function initialises the elementary mass matrix by taking
% advantage of the fields of the structure element given

function [M] = M_el_init(element,mat_prop)
% INPUTS
%   - element   : element to extract the matrix from        (structure)
%   - rho       : Material density                          (int)[kg.m^-3]

% OUTPUT
%   - M         : Elementary mass matrix of the element (cell array)

%% Extraction of the element properties
l = element.Length;         % Length                        (int)[m]
h = element.Height;         % Height                        (int)[mm]
b = element.Width;          % Width                         (int)[mm]
A = h*b*1e-6;               % Computing cross section area  (float)[m^2]

%% Extraction of the material properties
rho = mat_prop.rho;         % Density                       (int)[kg.m^-3]

%% Computation of the element properties needed
if h==b     
    % SQUARE BEAMS
    r = h/sqrt(12)*1e3;     % Computing radius of gyration  (float)[m]
else        
    % RECTANGULAR BEAM
    r = h/sqrt(12)*1e3;     % Computing radius of gyration  (float)[m]
    %% A VOIR PARCE QUE NORMALEMENT RADIUS OF GYRATION DEPENDS ON AXIS...
end

%% Initialisation of the elementary mass matrix

M = rho*A*l*...
[
1/3,         0,         0,     0,         0,         0, 1/6,         0,         0,      0,        0,         0; 
  0,     13/35,         0,     0,         0, -11*l/210,   0,      9/70,         0,      0,        0, -13*l/420;
  0,         0,     13/35,     0, -11*l/210,         0,   0,         0,      9/70,      0, 13*l/420,         0;
  0,         0,         0, r^2/3,         0,         0,   0,         0,         0, r^2/6,         0,         0; 
  0,         0, -11*l/210,     0,   l^2/105,         0,   0,         0, -13*l/420,      0, -l^2/140,         0;
  0,  11*l/210,         0,     0,         0,   l^2/105,   0,  13*l/420,         0,      0,        0,  -l^2/140;
1/6,         0,         0,     0,         0,         0, 1/3,         0,         0,      0,        0,         0;
  0,      9/70,         0,     0,         0,  13*l/420,   0,     13/35,         0,      0,        0, -11*l/210;
  0,         0,      9/70,     0, -13*l/420,         0,   0,         0,     13/35,      0, 11*l/210,         0; 
  0,         0,         0, r^2/6,         0,         0,   0,         0,         0,  r^2/3,        0,         0;
  0,         0,  13*l/420,     0,  -l^2/140,         0,   0,         0,  11*l/210,      0,  l^2/105,         0;
  0, -13*l/420,         0,     0,         0,  -l^2/140,   0, -11*l/210,         0,      0,        0,   l^2/105;
];
end