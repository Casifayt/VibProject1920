%% This function is used to plot a vector (mainly the fft)

% It plots the vector according to degrees of freedom of the
% accelerometer

% INPUTS : 
%   - q             : Vector to be plotted                          (array)[/]
%   - time_prop     : Time properties of the study                  (array)[s]
%   - elements_All  : Strucutre containing all elements properties  (struct)[/]
%   - nodes_All     : Matrix containing all nodes and their DOFs    (matrix)[/]

% OUTPUTS :
% /
function fourier_plot(q,time_prop, elements_All, nodes_All)

t_max = time_prop(1);
h = time_prop(2);

f_s = 1/h;
n_time_step = fix(t_max/h) + 1;
f_sampling = f_s * (0:(n_time_step/2))/n_time_step;

%% Extraction of accelerometer measurements
% Accelerometer is at node (4,15,3). This is the final node of the last element
% of the beam number 16, according to the indexation followed during the 
% discretisation (this is the beam going from (0,15,3) to (4,15,3)).

N_elem = numel(fieldnames(elements_All.('Beam1_elements')));

% Extracting the corresponding node number
accNode_nbr = elements_All.('Beam16_elements').(['Element' num2str(N_elem)]).nodeFin_nbr;

% Extracting the corresponding 3 translational DOFs from the nodes properties list
accDOFs_nbr = nodes_All(accNode_nbr,5:7);
i = 1;

% Calculating the dft of the displacement vector only for the DOFs of the
% accelerometer

for DOF = accDOFs_nbr(1):accDOFs_nbr(end)
    fft_q_disp{i} = fft(q(DOF,:), n_time_step);
    P2_disp{i} = abs(fft_q_disp{i}/n_time_step);
    P1_disp{i} = P2_disp{i}(1:ceil(n_time_step/2));
    P1_disp{i}(2:end-1) = 2*P1_disp{i}(2:end-1);
    i = i + 1;
end


    figure('Name',...
        ['Discrete Fourier Transform of q. Time properties t_max = ' num2str(t_max)...
        ' and h = ' num2str(h)]...
        ,'NumberTitle','off','Color','white'...
        ,'units','normalized','outerposition',[0 0 1 1]);
    
    title('Dynamic response according to the mode displacement method');
        
    subplot(3,1,1);
    plot(f_sampling, P1_disp{1}, 'color', [255,121,0]/256)
    title('X-dimension'); xlabel('Frequency [Hz]');ylabel('|FFT|');
    axis([0 10 0 Inf]); set(gca,'yscale','log');
%     [psor,lsor] = findpeaks(P1_disp{1}, f_sampling, 'SortStr','descend');
%     fprintf(' X : \n');
%     fprintf('Peak 1 : %.10f Hz\n Peak 2 : %.10f Hz\n Peak 3 : %.10f Hz\n',lsor(1), lsor(2),lsor(3));

    subplot(3,1,2);
    plot(f_sampling, P1_disp{2}, 'color', [255,121,0]/256)
    title('Y-dimension') ; xlabel('Frequency [Hz]');ylabel('|FFT|');
    axis([0 10 0 Inf]); set(gca,'yscale','log');
%     [psor,lsor] = findpeaks(P1_disp{2}, f_sampling, 'SortStr','descend');
%     fprintf(' Y : \n');
%     fprintf('Peak 1 : %.10f Hz\n Peak 2 : %.10f Hz\n Peak 3 : %.10f Hz\n',lsor(1), lsor(2),lsor(3));
    
    subplot(3,1,3);
    plot(f_sampling, P1_disp{3}, 'color', [255,121,0]/256)
    title('Z-dimension'); xlabel('Frequency [Hz]');ylabel('|FFT|');
    axis([0 10 0 Inf]); set(gca,'yscale','log');
%     [psor,lsor] = findpeaks(P1_disp{3}, f_sampling, 'SortStr','descend');
%     fprintf(' Z : \n');
%     fprintf('Peak 1 : %.10f Hz\n Peak 2 : %.10f Hz\n Peak 3 : %.10f Hz\n',lsor(1), lsor(2),lsor(3));

end