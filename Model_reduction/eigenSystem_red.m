%% Function resolving eigenvalue system and printing 8 first eigenmodes
% This function resolves the eigenvalue problem using the matlab function eig
% then extracts the eigenvalues which are the frequencies, in radians and
% then in hertz. It truncates those frequencies up to the tr-th one (tr for trunc).
% It may print the truncated frequencies. It also produces an array
% containg the corresponding rows in the matrix system for the sorted frequencies.

% INPUTS :
%   - K_red : Reduced stiffness matrix                              (matrix)[/]
%   - M_red : Reduced mass matrix                                   (matrix)[/]
%   - tr    : Number of frequencies desired                         (int)[/]
%   - f_fem : Natural frequencies computed by FEM method            (array)[Hz]

function eigenSystem_red(K_red, M_red, tr, f_fem)

% Initialisation of the timer
t_eigen = tic;

% Resolution of the eigenvalue system through MatLab eig function
[~, D] = eigs(K_red, M_red,tr+1, 'smallestabs');

% Extraction of exact frequencies in eigenvalues matrix D
wEx = sqrt((diag(D)));  % [rad/s]

% Sorting of the exact frequencies
wEx = sort(wEx);

% Initialising exact NX frequencies
NX = NX_exact_freq();

% Initialisation of counters for the print
i = 1;
eig_nbr = 1;
eig_freq_last = 0;
sum_err = 0;

% Choose between 
% f_fem to compare with frequencies from finite element method
% NX to compare with exact frequencies from NX
w_compare = f_fem;

fprintf('\nNatural frequencies after model reduction \n');
while eig_nbr < tr+1 && i < length(wEx)
   eig_freq_cur = wEx(i);
   if eig_freq_cur > eig_freq_last + .00001
       epsilon = abs(eig_freq_cur/2/pi - w_compare(eig_nbr))/w_compare(eig_nbr);
       sum_err = sum_err + epsilon;
       if 100 * epsilon < .01
            fprintf('%i : \x03C9 = %.3f rad;  f =  %.3f Hz, \x03B5 = %.2g %%\n',...
               eig_nbr, eig_freq_cur, eig_freq_cur/2/pi, 100 * epsilon);
       else
            fprintf('%i : \x03C9 = %.3f rad;  f =  %.3f Hz, \x03B5 = %.2f %%\n',...
               eig_nbr, eig_freq_cur, eig_freq_cur/2/pi, 100 * epsilon);
       end
       eig_freq_last = eig_freq_cur ;
       eig_nbr = eig_nbr + 1;
   end
   i = i + 1;
end

tElapsed = toc(t_eigen);
err_rel_avg = sum_err/tr;
if 100 * err_rel_avg < 1
    fprintf('Eigenvalue system solved in %.3fms, avg \x03B5 = %.2g %% \n',...
        1000 * tElapsed,100 * err_rel_avg);
else
    fprintf('Eigenvalue system solved in %.3fms, avg \x03B5 = %.2f %% \n',...
        1000 * tElapsed,100 * err_rel_avg);
end
end


function [NX_ex] = NX_exact_freq()
NX_ex = [
      3.88602597E-01;
      1.27103714E+00;
      2.17564722E+00;
      3.22022702E+00;
      3.80276261E+00;
      4.33157121E+00;
      4.34111495E+00;
      4.39099653E+00;
      4.49567716E+00;
      4.53934497E+00;
      4.56386850E+00;
      4.57658112E+00;
      4.79808884E+00;
      5.24877228E+00;
      6.28221899E+00;
      7.80863962E+00;
      8.22611650E+00;
      8.36882383E+00;
      8.48228180E+00;
      8.58073244E+00;
      8.76275621E+00;
      9.07867336E+00;
      9.21935507E+00;
      9.65851642E+00;
      9.84737019E+00;
      9.88552544E+00;
      1.01498932E+01;
      1.01572092E+01;
      1.01969226E+01;
      1.03375440E+01;
      1.04739154E+01;
      1.10088934E+01;
      1.10602803E+01;
      1.11360161E+01;
      1.11790002E+01;
      1.13686742E+01;
      1.18132508E+01;
      1.19570906E+01;
      1.19899536E+01;
      1.23159518E+01;
      1.25297097E+01;
      1.26200554E+01;
      1.26915677E+01;
      1.26986924E+01;
      1.27002511E+01;
      1.29023803E+01;
      1.29579876E+01;
      1.30649860E+01;
      1.33629407E+01;
      1.34834528E+01;
      ];
end