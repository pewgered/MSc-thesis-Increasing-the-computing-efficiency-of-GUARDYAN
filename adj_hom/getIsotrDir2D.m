function [v] = getIsotrDir2D()
x=sin(2*pi*rand);
v=[x,sqrt(1-x*x)];
end