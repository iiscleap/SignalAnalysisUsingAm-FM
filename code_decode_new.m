function C_res = code_decode_new(A,Fs,fp,flag_env_recon,flag_car_recon)
% ---------------------------------------------------------------------
% Usage C_res = code_decode_new(A,Fs,fp,flag_env_recon,flag_car_recon)
% ---------------------------------------------------------------------
% The function implements the FDLP decomposition and reconstruction.
% It uses a 32 band (bark-spaced) frequency decomposition and synthesis
% using QMF filters. The code also provides option to reconstruct the 
% signal directly from the envelope or the carrier signal alone.
% 
% Input options
% A              - signal in samples (req).
% Fs             - sampling frequency (8000).
% fp             - model order for FDLP analysis (20).
% flag_env_recon - flag to indicate if the output needs to be synthesized
%                  only based on the envelope. For speech signals, 
%                  this would result in whispered speech (0).
% flag_car_recon - flag to indicate if the output needs to be synthesized
%                  only based on the carrier. For speech signals, 
%                  this would result in contaminated speech without the message (0).
%
% Output
% C_res          - Resynthesized signal  
%
% Warning        - With default options output will be the same as input.
% -------------------------------------------------------------------------
% Sriram Ganapathy - Johns Hopkins University - 07-14-2011 - 
% See COPYING.TXT for GNU copyright and re-distribution.
% Some functions come from other colleagues. Acknowledgements to Petr
% Motlicek and Marios Athineos.
% -------------------------------------------------------------------------

if nargin < 1;  error('Oh Boy !!! Not Enough Input parameters');end
if nargin < 2 ; Fs = 8000; end
if nargin < 3 ; fp=20; end
if nargin < 4;  flag_env_recon = 0; end
if nargin < 5;  flag_car_recon = 0; end

scaling_factor = max(abs(A));         % To put the same gain to the encoded signal
siglen=length(A);

if scaling_factor < 1 
    flag_wav=1;
    A = A* 2^15;                      % wav to raw conversion.
    scaling_factor = scaling_factor*2^15;    
else
    flag_wav=0;
end
% ---------------------------------------------
% Basic Parameters
% ---------------------------------------------
sr = Fs;                             % This number found to have the minimum redundancy in MDCT calculation
nb = 32;                             % Number of QMF bands
flen = sr;                           % 1 sec window length (can be modified).
R = round(48*sr/1000);               % Overlap of frames to avoid blocking noise  

% -----------------------------------
% additional parameters:
% -----------------------------------

new_len = length(A) + (64 - mod(length(A),64));  % zero-padding of input sequence
A = [A; 2*randn(new_len - length(A),1)-1];       % RANDN ??again?? from (-1,1)

% this is needed for some examples when fft(1) coef. is too big and makes problems for next processing
A = A-mean(A); 
counter = 0;

% positions of selected input QMF outputs from X cell variable
p1 = [ones(1,26)*7 6 5 5 5 4 3];
p2 = [1:26, 14, 8, 9, 10, 6, 4];  

% --------------------------------------
% *********** Framing **************** %
% --------------------------------------
  

[fr_long, Nb_fr_long] = frame_new(A,flen,R);
wnd1 = [ones(1,flen-R), linspace(1,0,R)]';
wnd2 = [linspace(0,1,R) ones(1,flen-2*R), linspace(1,0,R)]';
wnd3 = [linspace(0,1,R) ones(1,flen-R)]';

for K = 1 :size(fr_long,2)
    if K == 1;
        fr_long(:,K) = fr_long(:,K).*wnd1;
    elseif K == Nb_fr_long
        fr_long(:,K) = fr_long(:,K).*wnd3;
    else
        fr_long(:,K) = fr_long(:,K).*wnd2;
    end
end


for K = 1:size(fr_long,2)
    disp(K); 
    x = fr_long(:,K);

% -----------------------------------
% sub-band decomposition
% -----------------------------------
[X] = QMF_7_analysis_W_bark_mod(x);
%   figure;
%   plot(X)
Y = cell(size(X));
 
% --------------------------------------------------------------------
% go over all bands and compute all pole model of temporal trajectory  
% --------------------------------------------------------------------
    
for I = 1:nb,
       
   counter = counter+1;   			     % for AA_q matrix    
   sig = X{p1(I),p2(I)};		     
    %N=5;
    N = round(fp*(flen/sr));
  
    SS = analysis_fdlp_burg(sig,N);     % Perform FDLP
    SS  = sqrt(SS'); 
    
    CARRIER = sig./SS;                  
    energy=sqrt(sum(sig.^2)/length(sig));
    
    % Envelope only reconstruction
    if flag_env_recon
        NEW_CARRIER = randn(size(CARRIER));
    else
        NEW_CARRIER = CARRIER;
    end
    
    % Carrier only reconstruction
    if flag_car_recon
       ENV = energy*ones(size(SS));
    else
        ENV = SS;
    end
    
    % Putting them together
    Y{p1(I),p2(I)} = NEW_CARRIER.*ENV;                % FDLP Synthesis !!!!!!!!!!

end

  B = QMF_7_synthesis_W_bark_mod(Y);                    % QMF synthesis from sub-bands
  if K == 1
     C_res = B';
  else
     ovr_region = C_res(end-R+1:end) + B(1:R)';        % Overlap add reconstruction
     C_res = [C_res(1:end-R);ovr_region; B(R+1:end)'];
  end
end;

C_res=scaling_factor*(C_res(1:siglen))/(max(abs(C_res(1:siglen))));    % 
if flag_wav 
    C_res = C_res/(2^15);
end

