%% Plotting function
% This function plots a 3D representation of the truss bridge with user
% inputs allowing to choose if we want to show the beams or elements nodes


function plotPerso(beams_All,elements_All)

% INPUTS
%   - beams_All     : Structure with all the beams      (structure)
%   - elements_All  : Structure with all the elements   (structure)

% OUTPUTS
%   - None (a graph)

% Extracting the number of elements
N_elem = numel(fieldnames(elements_All.Beam1_elements)); 

fct = input('\n Do you want a plot of the structure? Y/N\n','s');
if isempty(fct)
    return
else
    if fct == 'N' || fct == 'n'
        return
    else
        str = input('\n Do you want to plot the nodes ? Y/N\n','s');
        
        if str == 'Y' || str == 'y'
            str2 = input('\n Even the discretised elements nodes Y/N\n','s');
        end

        %% Plot of the total truss bridge
        for i = 1:numel(fieldnames(beams_All))
            beam = beams_All.(['Beam' num2str(i)]);
            P1 = beam.node_Initial;
            P2 = beam.node_Final;
                
            if str == 'Y' || str == 'y'
                plot3(P1(1),P1(2),P1(3),'o'); hold on;
                plot3(P2(1),P2(2),P2(3),'o');
            end
            pts = [P1;P2];                          
            plot3(pts(:,1), pts(:,2), pts(:,3));  grid on;
            if str ~= 'Y' && str ~= 'y'
                hold on;
            end
        end

        %% Plot of the nodes of the discretised elements
        
        for i = 1:numel(fieldnames(elements_All))
            current_beam = elements_All.(['Beam' num2str(i) '_elements']);
            for j = 1:numel(fieldnames(current_beam))
                current_element = current_beam.(['Element' num2str(j)]);
                P1 = current_element.node_Initial;
                P2 = current_element.node_Final;

                if str2 == 'Y' || str2 == 'y'
                    plot3(P1(1),P1(2),P1(3),'o');
                    plot3(P2(1),P2(2),P2(3),'o');
                end

                % Plot of each element seperately as a line
                %pts = [P1;P2];                       
                %plot3(pts(:,1), pts(:,2), pts(:,3));
            end
        end

        %% Plot settings
        view(54,27);            % sets user view accoridng to our view
        xlabel('x [m]');zlabel('z [m]'),ylabel('y [m]');
        title(['Truss bridge with each beam discretised in '  num2str(N_elem)  ' elements']);
    end
end
end