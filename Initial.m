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


%% Initiliasition of the beams
beams = beams_initialisation();

%% Plot of the truss bridge
for i = 1:length(beams)
    P1 = beams{1,i}{1,1};        
    P2 = beams{1,i}{1,2};        
    
    plot3(P1(1),P1(2),P1(3),'o');
    plot3(P2(1),P2(2),P2(3),'o');
    
    pts = [P1;P2];                          
    plot3(pts(:,1), pts(:,2), pts(:,3)); hold on; grid on;
    xlabel('[m]');zlabel('[m]'),ylabel('[m]');
end

%% Discretisation of the beams
N = 10;

for i = 1:length(beams)
    fprintf('\nDiscretisation of beam number %i into %i elements\n',i,N);
    elements = discretisation(beams{1,i},N,i);
    elementsAll{i,1}=elements;
end

%% Plot of the discretised elements
for i = 1:length(elementsAll)
    for j = 1:length(elementsAll{1,1})
        P1 = elementsAll{i,1}{1,j}{1,1};        
        P1 = elementsAll{i,1}{1,j}{1,2};

        plot3(P1(1),P1(2),P1(3),'o');
        plot3(P2(1),P2(2),P2(3),'o');

        pts = [P1;P2];                          
        %plot3(pts(:,1), pts(:,2), pts(:,3)); hold on; grid on;
    end
end
xlabel('[m]');zlabel('[m]'),ylabel('[m]');
title(['Truss bridge with each beam discretised in '  num2str(N)  ' elements']);



