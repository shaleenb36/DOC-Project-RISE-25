function [y] = curfit(m,wave)
% curve fit absorption data
		


abs=m(1)*exp(-(m(2)*wave))+m(3);

% abs=m(2)*wave+m(1);
y=abs;
