%% MECA0029 - Theory of Vibration
% Analysis of the dynamic behaviour of a truss bridge
% Academic year 2019 - 2020
%
% Alexandre ANDRIES and Casimir FAYT
% ULiege - Aerospatial Engineering

clear all;
format long;

%% Material properties
% Steel beams

E = 210e9;              % [Pa] Young's Modulus 
nu = .3;                % [1] Poissons's ratio
rho = 7.8*1e3;          % [kg.m^-3] Material density (steel, reference book p.390) 

%% Geometrial properties of the beams
% Lengths
sq3 = 3;                % [m] Square beams may be or 3 meters
sq4 = 4;                % [m] or 4 meters
sq5 = 5;                % [m] or 5 meters.
diagL = sqrt(sq5^2 + sq3^2); % [m] Diagonal beams length

% Diagonal beams are of rectangular section with 80 mm height x 30 mm width

diagH = 80*1e-3;        % [m] Diagonal beams height (rectangular section)
diagW = 30*1e-3;        % [m] Diagonal beams width (rectangular section)
areaDiag = diagH*diagW; % [m^2] Area of the rectangular beams

% Other beams are of square section of 70 mm edge
squareH = 70*1e-3;      % [m] Other beams edge (square section)
areaSquare = squareH^2; % [m^2] Area of the square beams


% Angle between beams
alphaRad = atan(3/5);   % [rad] Angle between horizontal and diagonal beams
alphaDeg = alphaRad*180/pi;
betaRad = atan(5/3);    % [rad] Angle between vertical and diagonal beams
betaDeg = betaRad*180/pi;

% Second moments of inertia
Isq3 = sq3^4/12; Isq4 = sq4^4/12; Isq5 = sq5^4/12;

Alexandre





