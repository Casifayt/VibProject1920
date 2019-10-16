%% MECA0029 - Theory of Vibration
% Analysis of the dynamic behaviour of a truss bridge
% Academic year 2019 - 2020
%
% Alexandre ANDRIES and Casimir FAYT
% ULiege - Aerospatial Engineering

clear all;
format short;

%% Material properties
% Steel beams

E = 210e9;              % [Pa] Young's Modulus 
nu = .3;                % [1] Poissons's ratio
rho = 7.8*1e3;          % [kg.m^-3] Material density (steel, reference book p.390) 


%% Initialisation of the beams
fprintf('\nInitialisation of the beams\n');
beams_All = beams_initialisation(0);
fprintf('\nBeams initialised\n');

%% Plot of the truss bridge
for i = 1:numel(fieldnames(beams_All))
    beam = beams_All.(['Beam' num2str(i)]);
    P1 = beam.node_Initial;
    P2 = beam.node_Final;
    
    plot3(P1(1),P1(2),P1(3),'o');
    plot3(P2(1),P2(2),P2(3),'o');
    
    pts = [P1;P2];                          
    plot3(pts(:,1), pts(:,2), pts(:,3)); hold on; grid on;
end
xlabel('x [m]');zlabel('z [m]'),ylabel('y [m]');
clear beam;

%% Discretisation of the beams
N = 10;
fprintf('\nDiscretisation of the beams into %i elements\n',N);
for i = 1:numel(fieldnames(beams_All))
    beam = beams_All.(['Beam' num2str(i)]);
    elements = discretisation(beam,N,i,0);
    elements_All.(['Beam' num2str(i) '_elements'])=elements;
end
fprintf('\nBeams discretised\n');


%% Plot of the nodes of the discretised elements
for i = 1:numel(fieldnames(elements_All))
    current_beam = elements_All.(['Beam' num2str(i) '_elements']);
    for j = 1:numel(fieldnames(current_beam))
        current_element = current_beam.(['Element' num2str(j)]);
        P1 = current_element.node_Initial;
        P2 = current_element.node_Final;

        plot3(P1(1),P1(2),P1(3),'o');
        plot3(P2(1),P2(2),P2(3),'o');

        pts = [P1;P2];                          
        %plot3(pts(:,1), pts(:,2), pts(:,3)); hold on; grid on;
    end
end
title(['Truss bridge with each beam discretised in '  num2str(N)  ' elements']);

%% Cleansing of the useless variables left
var_2_clean = {'current_beam', 'current_element', 'i', 'j','N','beam', 'P1', 'P2', 'pts', 'var_2_clean'};
clear (var_2_clean{:});