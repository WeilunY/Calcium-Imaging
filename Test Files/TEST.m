function [  ] = TEST(  )
%TEST Summary of this function goes here
%   Detailed explanation goes here

x = linspace(0,2,1100);

a = 100;
b = 1;
y = 1./(1+a*exp(-x)).^b;

plot(x,y)

end

