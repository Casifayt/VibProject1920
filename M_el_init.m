%% Function initialising the elementary mass matrix
% This function initialises the elementary mass matrix by taking
% advantage of the fields of the structure element given

function [M] = M_el_init(element,rho)
% INPUTS
%   - element   : element to extract the matrix from        (structure)
%   - rho       : Material density                          (int)[kg.m^-3]

% OUTPUT
%   - M         : Elementary mass matrix of the element (cell array)

%% Extraction of the element properties
l = element.Length;         % Extracting element length (l)
h = element.Height;         % Extracting element height (mm)
b = element.Width;          % Extracting element width (mm)
A = h*b;                    % Computing cross section are (mm^2)

%% Computation of the element properties needed
if h==b     
    % SQUARE BEAMS
    r = h/sqrt(12)*1e3;     % Computing radius of gyration (m)
else        
    % RECTANGULAR BEAM
    r = h/sqrt(12)*1e3;     % Computing radius of gyration (m)
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