function [B] = QMF_7_synthesis_W_bark_mod(X)
warning off;

% Thanks to Petr Motlicek, Idiap Research Institute for the QMF filter
% design.

% X - matrix with signal decomposition into 5 levels == 32 branches
% X{1,1} = A;
% X{2,... X{7,64}

% B - synthesized signal
% Y - model variable
% Y{7,64} ... Y{1,1}
% Y{1,1} - signal

load mfilters2.mat ;
F = cell(1,2);

Y = [];

Y = cell(7,2^6);

% First level
for i = 1:26,
  Y{7,i} = zeros(length(X{7,i})*2,1);
  Y{7,i}(1:2:end) = X{7,i};
end;
N = 99;
% Second level
F{1} = GlpF6;
F{2} = GhpF6;
N6 = 13;
k = 1:2:64;
for i = 1:13,  
  Y{6,i} = circonv_mod(Y{7,k(i)},F{1},length(Y{7,k(i)})) + circonv_mod(Y{7,k(i)+1},F{2},length(Y{7,k(i)+1}));
  Y{6,i} = [Y{6,i}(N6+1:end) Y{6,i}(1:N6)];             % Circular shift due to circular convolution
  Hlp = Y{6,i};
  Y{6,i} = zeros(length(Y{6,i})*2,1);
  Y{6,i}(1:2:end) = Hlp;
end;


Y{6,14} = zeros(length(X{6,14})*2,1); 
%size(Y{6,14}(1+99*2+1*2:2:end))
%size(X{6,14}(1+99+1:end))
Y{6,14}(1:2:end) = X{6,14}(1:end);
%Y{6,14}(1+(99)*2:2:end) = X{6,14}(1:end-(99));
  
 
% Third level 

F{1} = GlpF5;
F{2} = GhpF5;
N5 = 17;

k = 1:2:32;
for i = 1:7, 
  Y{5,i} = circonv_mod(Y{6,k(i)},F{1},length(Y{6,k(i)})) + circonv_mod(Y{6,k(i)+1},F{2},length(Y{6,k(i)+1}));
  Y{5,i} = [Y{5,i}(N5+1:end) Y{5,i}(1:N5)];             % Circular shift due to circular convolution  
  Hlp = Y{5,i};
  Y{5,i} = zeros(length(Y{5,i})*2,1);
  Y{5,i}(1:2:end) = Hlp;
end;

Y{5,8} = zeros(length(X{5,8})*2,1); 
Y{5,8}(1:2:end) = X{5,8};
%Y{5,8}(1+(297)*2:2:end) = X{5,8}(1:end-297);

Y{5,9} = zeros(length(X{5,9})*2,1); 
Y{5,9}(1:2:end) = X{5,9};
%Y{5,9}(1+(297)*2:2:end) = X{5,9}(1:end-297);

Y{5,10} = zeros(length(X{5,10})*2,1); 
Y{5,10}(1:2:end) = X{5,10};
%Y{5,10}(1+(297)*2:2:end) = X{5,10}(1:end-297);

% Fourth level
F{1} = GlpF4;
F{2} = GhpF4;
N4 = 25;

k = 1:2:16;
for i = 1:5,  
  Y{4,i} = circonv_mod(Y{5,k(i)},F{1},length(Y{5,k(i)})) + circonv_mod(Y{5,k(i)+1},F{2},length(Y{5,k(i)+1}));
  Y{4,i} = [Y{4,i}(N4+1:end) Y{4,i}(1:N4)];             % Circular shift due to circular convolution    
  Hlp = Y{4,i};
  Y{4,i} = zeros(length(Y{4,i})*2,1);
  Y{4,i}(1:2:end) = Hlp;
end;

Y{4,6} = zeros(length(X{4,6})*2,1); 
Y{4,6}(1:2:end) = X{4,6};
%Y{4,6}(1+(693)*2:2:end) = X{4,6}(1:end-693);


% Fifth level
F{1} = GlpF3;
F{2} = GhpF3;
N3 = 51;

k = 1:2:8;
for i = 1:3,  
  Y{3,i} = circonv_mod(Y{4,k(i)},F{1},length(Y{4,k(i)})) + circonv_mod(Y{4,k(i)+1},F{2},length(Y{4,k(i)+1}));
  Y{3,i} = [Y{3,i}(N3+1:end) Y{3,i}(1:N3)];             % Circular shift due to circular convolution
  Hlp = Y{3,i};
  Y{3,i} = zeros(length(Y{3,i})*2,1);
  Y{3,i}(1:2:end) = Hlp;
end;

Y{3,4} = zeros(length(X{3,4})*2,1); 
Y{3,4}(1:2:end) = X{3,4};
%Y{3,4}(1+(1485)*2:2:end) = X{3,4}(1:end-1485);


% Sixth level
load mfilters.mat ;
F{1} = GlpF;
F{2} = GhpF;
N2 = 99;
k = 1:2:4;
for i = 1:2,  
  Y{2,i} = circonv_mod(Y{3,k(i)},F{1},length(Y{3,k(i)})) + circonv_mod(Y{3,k(i)+1},F{2},length(Y{3,k(i)+1}));
  Y{2,i} = [Y{2,i}(N2+1:end) Y{2,i}(1:N2)];             % Circular shift due to circular convolution  
  Hlp = Y{2,i};
  Y{2,i} = zeros(length(Y{2,i})*2,1);
  Y{2,i}(1:2:end) = Hlp;
end;

F{1} = GlpF;
F{2} = GhpF;
N1 = 99;

Y{1,1} = circonv_mod(Y{2,1},F{1},length(Y{2,1})) + circonv_mod(Y{2,2},F{2},length(Y{2,2}));
Y{1,1} = [Y{1,1}(N1+1:end) Y{1,1}(1:N1)];             % Circular shift due to circular convolution
% shift 1455 samples for 48kHz == 3.03ms
%B = Y{1,1}(1455:end);
%B = Y{1,1}(3070:end);
B = Y{1,1};
