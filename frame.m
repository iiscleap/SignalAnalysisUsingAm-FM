function [xtr,NB_FR]=frame(x,FRAME,R)
% description : [xtr,NB_FR]=frame(x,FRAME,R)
% sructuration of signal x to the  NB_FR of frames
% FRAME is length of frames which you want 
% R is overlapping of frames
NB_FR=fix((length(x)-FRAME)/(FRAME-R)+1);
[l,c]=size(x);
if(l>c)
x=x';
end

if (NB_FR ~= 0)
for tr=1:NB_FR
	xtr(:,tr)= x(1+(tr-1)*(FRAME-R):tr*(FRAME-R)+R)';
end
else disp(' can not create frames')
end


