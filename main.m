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

E = 210e9;              % Young's Modulus       [Pa] 
nu = .3;                % Poissons's ratio      [/]
G = .5*E/(1+2*nu);      % Shear modulus         [Pa]  
rho = 7.8*1e3;          % Material density      [kg.m^-3] 
% Density for steel extracted from reference book p.390

mat_prop = struct('E',E,'nu',nu,'rho',rho,'G',G);
% All material properties stored in a structure with corresponding fields


%% Initialisation of the beams
fprintf('\nInitialisation of the beams\n');
print = 0;              % If 1, prints details beams_initialisation
                        % no print if 0
beams_All = beams_initialisation(print);
fprintf('\nBeams initialised\n');

%% Discretisation of the beams
N_elem = 10;         % Number of elements per beam
fprintf('\nDiscretisation of the beams into %i elements\n',N_elem);
for i = 1:numel(fieldnames(beams_All))
    beam = beams_All.(['Beam' num2str(i)]);
    elements = discretisation(beam,N_elem,i,0);
    elements_All.(['Beam' num2str(i) '_elements'])=elements;
end
fprintf('\nBeams discretised\n');


%% Plotting according to user inputs
plotPerso(beams_All,elements_All);

%% Initialisation of elementary matrices for each element
for i = 1:numel(fieldnames(elements_All))           % Scanning on each beam
    current_beam = elements_All.(['Beam' num2str(i) '_elements']);  % Selecting one beam
    
    for j = 1:numel(fieldnames(current_beam))       % Scanning on each element of this beam
        current_element = current_beam.(['Element' num2str(j)]);    % Selecting one element
        K_el = K_el_init(current_element,mat_prop);      % Initialising the elementary stiffness matrix
        M_el = M_el_init(current_element,mat_prop);      % Initialising the elementary mass matrix
    end
end



%% Localisation matrix (L. III sl.17)
locel = zeros(N_elem,12);





%% Cleansing of the useless variables left
var_2_clean = {'current_beam', 'nodes_beam_plotting','nodes_elements_ploting',...
    'current_element', 'print','i', 'j','N','beam', 'P1', 'P2', 'pts', ...
    'elements','var_2_clean'};
clear (var_2_clean{:});