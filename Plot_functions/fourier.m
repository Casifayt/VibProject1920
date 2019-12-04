%% Function calculating the discrete fast fourier transform of a vector
% This function computes the discrete fft by using the MatLab built-in fft
% function. The time properties of the problem are used to define sampling
% frequency.

% INPUTS :
%   - q         : Vector to undergo fft                     (array)[/]
%   - time_prop : Array containing time properties          (array)[s]

% OUTPUTS :
%   - q_fft     : Discrete fft of the vector q              (array)[/]
%   - f_sampling: Sampling frequency corresponding to fft   (array)[Hz]

function [q_fft, f_sampling] = fourier(q,time_prop)

% Extracting time properties
t_max = time_prop(1);
h = time_prop(2);

% Computing sampling frequency
f_sampling = 1/h;

% Computing total number of discrete time steps
n_time_step = t_max / h;

% Computing the fft
q_fft = fft(q,n_time_step);

% Restricting to half the spectrum (other half not of interest)
q_fft = q_fft(1:n_time_step/2);

% Computing sampling frequency according to spectrum
f_sampling = f_sampling*(0:length(q_fft) - 1)/n_time_step;

end