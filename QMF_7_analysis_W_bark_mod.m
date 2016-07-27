function [X] = QMF_7_analysis_W_bark_mod(A)
warning off;
% Function to perform the QMF analysis with low delay and with sharp
% transistion bandwidth
% Output contains a cell array with size (7,64) with outputs from each
% stage in a 6 stage binary decomposition into 64 bands. Appropiate bands
% can be chosen to obtain 32 bands which are non-uniform.
% Sriram Ganapathy - Johns Hopkins University - 07-14-2011.

% Thanks to Petr Motlicek, Idiap Research Institute for the QMF filter
% design.

% parameter setting
load mfilters.mat;

%load /home/speech/pmotlic/AMR_WB+_examples/FIR_LP;

H = cell(1,2);
F = H;
H{1} = HlpF ;
H{2} = HhpF;

% 
X = [];
X = cell(7,2^6);
X{1,1} = A';

% First level
% sampling at 24kHz
for i = 1:2,
  X{2,i} = circonv_mod(X{1,1},H{i},length(X{1,1})); 
  X{2,i} = X{2,i}(1:2:end);
end;

% Second level
% sampling at 12kHz
H{1} = HlpF;
H{2} = HhpF;

k = [1 1 2 2];
for i = 1:4,
  X{3,i} = circonv_mod(X{2,k(i)},H{mod(i-1,2)+1},length(X{2,k(i)}));
  X{3,i} = X{3,i}(1:2:end);
end;
%X{3,4} = B32
% ---------------------------------------------------------
% Third level
% sampling at 6kHz
load mfilters2.mat;
H{1} = HlpF3;
H{2} = HhpF3;

k = [1 1 2 2 3 3 4 4];
for i = 1:8,
  X{4,i} = circonv_mod(X{3,k(i)},H{mod(i-1,2)+1},length(X{3,k(i)}));
  X{4,i} = X{4,i}(1:2:end);
end;
%X{4,6} = B31

% Fourth level
% sampling at 3kHz
H{1} = HlpF4;
H{2} = HhpF4;

k = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
for i = 1:16,
  X{5,i} = circonv_mod(X{4,k(i)}, H{mod(i-1,2)+1},length(X{4,k(i)}));
  X{5,i} = X{5,i}(1:2:end);
end;
% Fifth level
% sampling at 1.5kHz
H{1} = HlpF5;
H{2} = HhpF5;

k = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16];
for i = 1:32,
  X{6,i} = circonv_mod(X{5,k(i)},H{mod(i-1,2)+1},length(X{5,k(i)}));
  X{6,i} = X{6,i}(1:2:end);
end;

% ---------------------------------------------------------
% Sixth level
% sampling at 750Hz
H{1} = HlpF6;
H{2} = HhpF6;

k = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 ...
17 17 18 18 19 19 20 20 21 21 22 22 23 23 24 24 25 25 26 26 27 27 28 28 29 29 30 30 31 31 32 32];
for i = 1:64,
  X{7,i} = circonv_mod(X{6,k(i)},H{mod(i-1,2)+1},length(X{6,k(i)}));
  X{7,i} = X{7,i}(1:2:end);
end;
%X{7,1} = B1
%...
%X{7,26} = B26
