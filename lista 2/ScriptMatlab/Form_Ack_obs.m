function [ L ] = Form_Ack_obs(A,C,u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
syms s
Pol=(s-u(1));
for i=2:length(u)
    Pol=Pol*(s-u(i));
end
[f c]=size(A);
M1=zeros(f,1);
M1(end)=1;
Pol=simplify(Pol);
Pol=collect(Pol);
Vec1=sym2poly(Pol);
Vec1=fliplr(Vec1);
Pobs_A=zeros(size(A));
for i=1:length(Vec1)
    Pobs_A=Pobs_A+Vec1(i)*A^(i-1);
end
Ok=[C];
for i=1:f-1
    Ok=[Ok;C*A^i];
end
L=Pobs_A*inv(Ok)*M1;

end

