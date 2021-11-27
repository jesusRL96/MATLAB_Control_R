function [ K ] = Funcion_FormAckerman( A,B,u )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
syms s
Pol=(s-u(1));
for i=2:length(u)
    Pol=Pol*(s-u(i));
end
[f c]=size(A);
M1=[1:f]*0;
M1(end)=1;
Pol=simplify(Pol);
Pol=collect(Pol);
Vec1=sym2poly(Pol);
Vec1=fliplr(Vec1);
PLC_A=zeros(size(A));
for i=1:length(Vec1)
    PLC_A=PLC_A+Vec1(i)*A^(i-1);
end
CK=[B];
for i=1:f-1
    CK=[CK,A^i*B];
end
K=M1*inv(CK)*PLC_A;
end

