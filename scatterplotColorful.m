function [ ] = scatterplotColorful( y,x,yeq,titleStr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
yvec = y(:);
xvec = x(:);
yzfvec = yeq(:);

hzf = scatterplot(yvec(abs(xvec-(-1))<1E-4),1,0,'or'); hold on;
scatterplot(yvec(abs(xvec-(1))<1E-4),1,0,'xr',hzf); 
scatterplot(yzfvec(abs(xvec-(-1))<1E-4),1,0,'ob',hzf);
scatterplot(yzfvec(abs(xvec-(1))<1E-4),1,0,'xb',hzf);
hold off;
title(titleStr);

end

