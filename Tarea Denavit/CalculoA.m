function [ A ] = CalculoA(D)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
a=D(1);     alpha=D(2);         d=D(3);     theta=D(4);
a=sym(a);   alpha=sym(alpha);   d=sym(d);   theta=sym(theta);
I=sym([a,alpha,d,theta]);
A=[cos(theta),-sin(theta)*cos(alpha),sin(theta)*sin(alpha),a*cos(theta);...
    sin(theta),cos(theta)*cos(alpha),-cos(theta)*sin(alpha),a*sin(theta);...
    0,sin(alpha),cos(alpha),d;...
    0,0,0,1];
A=simplify(A);
end

