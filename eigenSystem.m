%% Function resolving eigenvalue system and printing 8 first eigenmodes
% This function resolves the eigenvalue problem using the matlab function eig
% then extracts the eigenvalues which are the frequencies, in radians and
% then in hertz. It truncates those frequencies up to the tr-th one (tr for trunc).
% It may print the truncated frequencies. It also produces an array
% containg the corresponding rows in the matrix system for the sorted frequencies.

% INPUTS :
%   - K_S : Structural stiffness matrix                             (matrix)[/]
%   - M_S : Structural mass matrix                                  (matrix)[/]
%   - tr  : Number of frequencies desired                           (int)[/]

% OUTPUTS :
%   - fEx   : Sorted exact eigen frequencies                        (array)[Hz]
%   - f_tr  : Sorted truncated up to the trunc-th one eigen freq    (array)[Hz]
%   - wEx   : Sorted exact eigen angular frequencies                (array)[rad]
%   - w_tr  : Sorted truncated up to the trunc-th one eigen ang freq(array)[rad]
%   - wColInV : Corresponding row of the sorted frequencies in D    (array)[/]
%   - V     : Full matrix solution of the eigenvalue problem        (matrix)[/]
%   - D     : Diagonal matrix solution of the eigenvalue problem    (matrix)[/]

function [fEx, f_tr,wEx,w_tr,wColInV,V,D] = eigenSystem(K_S,M_S,tr)

% Initialisation of the timer
t_eigen = tic;

w_tr = zeros(1,tr);

% Resolution of the eigenvalue system through MatLab eig function
[V, D]=eig(K_S, M_S);

% Extraction of exact frequencies in eigenvalues matrix D
wEx = sqrt(diag(D));  % [rad/s]

% Sorting of the exact frequencies
wEx = sort(wEx);

% Initialisation of the exact frequencies produced with Siemens NX
NX_ex = NX_exact_freq();

% Initialisation of counters for the print
i = 1;
eig_nbr = 1;
eig_freq_last = 0;

fprintf(['\nNatural frequencies of the ' num2str(tr) ' first eigenmodes : \n']);
while eig_nbr < tr+1 && i < length(wEx)
   eig_freq_cur = wEx(i);
   if eig_freq_cur > eig_freq_last + .00001
       epsilon = (eig_freq_cur - NX_ex(eig_nbr))/NX_ex(eig_nbr);
       fprintf('%i : \x03C9 = %f rad;  f =  %f Hz, \x03B5 = %.2f %%\n',...
           eig_nbr, eig_freq_cur, eig_freq_cur/2/pi,epsilon);
       eig_freq_last = eig_freq_cur ;
       w_tr(eig_nbr) = eig_freq_cur;
       eig_nbr = eig_nbr + 1;
   end
   i = i + 1;
end

tElapsed = toc(t_eigen);
fprintf('Eigenvalue system solved in %.2fs\n',tElapsed);


fEx = wEx /2/pi;
f_tr = w_tr/2/pi;   

w_temp = sqrt(diag(D));
wColInV = zeros(1,length(w_temp));

for i = 1:length(w_temp)
    for j = 1:length(wEx)
       if w_temp(i) ==  wEx(j)
           if wColInV(end,j) == 0
               wColInV(end,j) = i;
           else
               wColInV(end+1,j) = i;
           end
       end
   end
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