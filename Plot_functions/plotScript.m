%% Plotting function
% This function plots a 3D representation of the truss bridge with user
% inputs allowing to choose if we want to show the beams or elements nodes

% INPUTS
%   - beams_All     : Structure with all the beams      (structure)
%   - elements_All  : Structure with all the elements   (structure)
%   - tr            : Number of truncated frequencies   (int)[/]
%   - wColInV : Corresponding row of the sorted frequencies in D    (array)[/]
%   - V             : Structural displacements matrix   (cell array)[/]
%   - locel         : Localisation element matrix       (cell array)[/]

% OUTPUTS : / (a plot)


function plotScript(beams_All,elements_All, tr, wColInV, V, locel)

% Extracting the number of elements
N_elem = numel(fieldnames(elements_All.Beam1_elements)); 

    %% Plot of the deformed structure

    % Choice of the mode
    mode = input(['\nWhat mode (smaller than ' num2str(tr+1) ')?\n']);

    % Error management
    if mode > tr
        error(['I asked you to choose a mode smaller than ' num2str(tr+1) '!']);
    end

    % Construction of deformed elements
    elements_Deformed = deformation(mode,wColInV,V,locel,elements_All);

    % If one wants to plot the initial discretised structure
    % elements_Deformed = elements_All;

    for i = 1:numel(fieldnames(elements_Deformed))

        current_beam = elements_Deformed.(['Beam' num2str(i) '_elements']);

        for j = 1:numel(fieldnames(current_beam))

            current_element = current_beam.(['Element' num2str(j)]);

            P1 = current_element.node_Initial;
            P2 = current_element.node_Final;
            pts = [
                P1;
                P2;
                ];
            plot3(pts(:,1), pts(:,2), pts(:,3),'Color','r'); hold on; grid on;
        end
    end

    %% Plot of the undeformed truss bridge
    for i = 1:numel(fieldnames(beams_All))

        beam = beams_All.(['Beam' num2str(i)]);

        P1 = beam.node_Initial;
        P2 = beam.node_Final;
        pts = [
            P1;
            P2;
            ];

        plot3(pts(:,1), pts(:,2), pts(:,3),'--','Color','k');
    end


    %% Plot settings
    view(54,27);            % sets user view accoridng to our view
    xlabel('x [m]');zlabel('z [m]'),ylabel('y [m]');
%     axis([-1 5 -1 21 -1 4]);
    axis equal;
    if mode == 1 
        str = 'st';
    elseif mode == 2
        str = 'nd';
    elseif mode == 3
        str = 'rd';
    else 
        str = 'th';
    end
%     title({'Truss bridge', ['Discretisation in '  num2str(N_elem)...
%         ' elements, ' num2str(mode) str ' mode']});
    set(gca,'ColorScale','log') 
end