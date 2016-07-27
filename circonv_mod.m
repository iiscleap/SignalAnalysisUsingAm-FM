%CIRCONV	N-point circular convolution
%
%	C = CIRCONV(A,B,N) performs the N-point circular convolution
%	of vectors A and B.  C is returned as a row vector.  A and B
% 	must be vectors, but may be of different lengths.  N must be
% 	a positive, non-zero integer.  The results of CIRCONV will
%	match that of CONV if N>=( length(A) + length(B) - 1).  This
%	is also the alias-free criterion for the circular convolution.
%
%	See also CONV

%
% Edward Brian Welch
% edwardbrianwelch@yahoo.com
% 16-Feb-1999
%
function[C] = circonv_mod(A,B,N)

if size(A,1) ~=1
    A = transpose(A);
end
if size(B,1) ~=1
    B = transpose(B);
end


C = real(ifft(fft(A,N).*fft(B,N)));
