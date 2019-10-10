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

E = 210*1e9;            % [Pa] Young's Modulus 
nu = .3;                % [1] Poissons's ratio
rho = 8*1e3;            % [kg.m^-3] Material density (steel)

%% Geometrial properties of the beams
% Diagonal beams are of rectangular section with 80 mm height x 30 mm width

diagH = 80*1e-3;        % [m] Diagonal beams height (rectangular section)
diagW = 30*1e-3;        % [m] Diagonal beams width (rectangular section)

% Other beams are of square section of 70 mm edge
beamSq = 70*1e-3;       % [m] Other beams edge (square section)

% Angle between beams
alphaRad = atan(3/5);   % [rad] Angle between horizontal and diagonal beams
alphaDeg = alphaRad*180/pi;
betaRad = atan(5/3);    % [rad] Angle between vertical and diagonal beams
betaDeg = betaRad*180/pi;






