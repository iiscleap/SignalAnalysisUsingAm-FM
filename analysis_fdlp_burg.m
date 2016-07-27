function hilb_env_new = analysis_fdlp_burg(x,p)

% To perform FDLP to analyze the sub-band signal.
% compute to AR model of the Hilbert Envelope of the signal.
% The derivation is explained in S. Ganapathy, et al, "Autoregressive 
% Models of Amplitude Modulations In Audio Compression",IEEE TASLP, 2010.
% Sriram Ganapathy - Johns Hopkins University - 07-14-2011. 

if size(x,2) == 1
    x = x';                                     % Make a column vector
end

N = length(x);
M = 2*N - 1;                                    % Size of the even-symmetrized signal.

x_e = [x fliplr(x(2:end))];                     % Even component of the signal.
x_fft = real(fft(x_e,M));                       % Real FFT
x_fft = x_fft(1:N);                             % Keeping N real samples of the FFT
% x_fft = x_fft.*hanning(N)';                   % Uncomment for freq. window.

y = [x_fft(1) 2*x_fft(2:end)];                  % Analytic signal - one side FFT

%-------------------------------------------------------------------
%        Auto-regressive modelling of Hilbert Envelopes
% --------------------------------------------------------------------

hilb_env_new = pyulear(y/sqrt(N),p,M);

 


