%% Function initialising the elementary mass matrix
% This function initialises the elementary mass matrix by taking
% advantage of the fields of the structure element given

function [M] = M_el_init(element,mat_prop)
% INPUTS
%   - element   : element to extract the matrix from        (structure)
%   - rho       : Material density                          (int)[kg.m^-3]

% OUTPUT
%   - M         : Elementary mass matrix of the element (cell array)

%% Extraction of the element geometric properties
l = element.Length;         % Length                        (int)[m]
h = element.Height*1e-3;    % Height                        (int)[m]
b = element.Width*1e-3;     % Width                         (int)[m]
A = h*b;                    % Computing cross section area  (float)[m^2]

%% Extraction of the element material properties
rho = mat_prop.rho;         % Density                       (int)[kg.m^-3]

%% Computation of the element geometric properties needed
if h==b     
    % SQUARE BEAMS
    I = h^4/12;             % Second moment of area         (float)[m^4]
    r = sqrt(I/A);          % Computing radius of gyration  (float)[m]
else        
    % RECTANGULAR BEAM
    Iy = h*b^3/12;          % Smallest second moment of area(float)[m^4]
    % Only the smallest one is taken because it is the one corresponding to
    % the smallest radius of gyration, which corresponds to the  likeliest 
    % axis to suffer from buckling, giving therefore the limit value
    
    r = sqrt(Iy/A);         % Computing radius of gyration  (float)[m]
end

%% Computation of the elementary mass matrix

M = rho*A*l*...
[
1/3,         0,         0,     0,         0,         0, 1/6,         0,         0,      0,        0,         0; 
  0,     13/35,         0,     0,         0, -11*l/210,   0,      9/70,         0,      0,        0, -13*l/420;
  0,         0,     13/35,     0, -11*l/210,         0,   0,         0,      9/70,      0, 13*l/420,         0;
  0,         0,         0, r^2/3,         0,         0,   0,         0,         0,  r^2/6,        0,         0; 
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