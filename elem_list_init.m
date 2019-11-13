function [elem_list] = elem_list_init(elements_All)
    
%elem_list = zeros(N_elem,2);

num_beams = numel(fieldnames(elements_All)); 
num_el = numel(fieldnames(elements_All.Beam1_elements));
unit_vec_x = [1;0;0];

for i = 1: num_beams                % Scanning on each beam
    current_beam = elements_All.(['Beam' num2str(i) '_elements']);  % Selecting one beam
    if num_el > 1
        for j =1:num_el
            current_element = current_beam.(['Element' num2str(j)]);    % Selecting one element
            if current_element.Orientation*unit_vec_x == 1
               %elem_list(nbr_el_done,1) = 
            end
        end        
    end 
end
end