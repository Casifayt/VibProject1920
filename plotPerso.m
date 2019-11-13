%% Plotting function
% This function plots a 3D representation of the truss bridge with user
% inputs allowing to choose if we want to show the beams or elements nodes

% INPUTS
%   - beams_All     : Structure with all the beams      (structure)
%   - elements_All  : Structure with all the elements   (structure)
%   - tr            : Number of truncated frequencies   (int)[/]

% OUTPUTS : / (a plot)


function plotPerso(beams_All,elements_All, tr)

% Extracting the number of elements
N_elem = numel(fieldnames(elements_All.Beam1_elements)); 

fct = input('\nPlot of the deformed structure? Y/N\n','s');
if isempty(fct)
    return
else
    if strcmp(fct, 'N') == 1 || strcmp(fct, 'n') == 1
        return
    else
        
        mode = input(['What mode (smaller than ' num2str(tr+1) ')?\n']);
        
        if mode > tr
            error(['I asked you to choose a mode smaller than ' num2str(tr+1) '!']);
            return;
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
            
            plot3(pts(:,1), pts(:,2), pts(:,3),'--','Color','k');  hold on; grid on;
        end

        %% Plot of the discretised deformed elements
        
        for i = 1:numel(fieldnames(elements_All))
            
            current_beam = elements_All.(['Beam' num2str(i) '_elements']);
            
            for j = 1:numel(fieldnames(current_beam))
                
                current_element = current_beam.(['Element' num2str(j)]);
                
                P1 = current_element.node_Initial;
                P2 = current_element.node_Final;
                pts = [
                    P1;
                    P2;
                    ];
                
                plot3(pts(:,1), pts(:,2), pts(:,3),'Color','r');                
            end
        end

        %% Plot settings
        view(54,27);            % sets user view accoridng to our view
        xlabel('x [m]');zlabel('z [m]'),ylabel('y [m]');
        axis([-1 5 -1 21 -1 4]);
        if mode == 1 
            str = 'st';
        elseif mode == 2
            str = 'nd';
        elseif mode == 3
            str = 'rd';
        else 
            str = 'th';
        end
        title({'Truss bridge', ['Discretisation in '  num2str(N_elem)...
            ' elements, ' num2str(mode) str ' mode']});
    end
end
end