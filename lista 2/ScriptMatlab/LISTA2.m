clc
clear all
close all
syms k Bl Bm Jl Jm theta_l theta_m
A=[0 1 0 0;-k/Jl -Bl/Jl k/Jl 0;0 0 0 1;k/Jm 0 -k/Jm -Bm/Jm];
B=[0;0;0;1/Jm];
C=[1 0 0 0];
syms s
% Ejercicio 1.2 (Ganancia)<<Gs=C(sI-A)B>>
disp('Ejercicio 1.2 (Ganancia)<<Gs=C(sI-A)B>>')
Gs=C*inv(eye(size(A))*s - A)*B

%% Ejercicio 2
disp('EJERCICIO 2')
k=.8; Bl=.01; Bm=.015; Jl=.0004; Jm=.0004;      % Valores numericos
syms z1 z2 wn1 wn2
P_d=(s^2+2*z1*wn1*s+wn1^2)*(s^2+2*z2*wn2*s+wn2^2);   % Polinomio de dinamica deseada
z1=1; z2=2; wn1=10; wn2=15;     % Valores numericos
P_d=subs(P_d);
disp('Valores propios de la dinamica deseada:')
u=solve(P_d)      % Valores propios de la dinamica deseada
A=subs(A);
B=subs(B);
disp('Constantes de retroalimentacion:')
K=double(Funcion_FormAckerman(A,B,u))   % Calculo de constantes de retroalimentacion
disp('Comprobacion (comp=eig(A-B*K) ):')
comp=eig(A-B*K)     % Comprobacion 
%% Ejercicio 3
disp('EJERCICIO 3')
disp('Constantes para observabilidad:')
syms z1o z2o wn1o wn2o
u_o=[-40 -40 -40 -40];          % Polos propuestos para observabilidad (a la izquierda de los polos de dinamica deseada)
P_do=poly(u_o)*[s^4;s^3;s^2;s;1];
P_de=(s^2+2*z1o*wn1o*s+wn1o^2)*(s^2+2*z2o*wn2o*s+wn2o^2);
P_doC=coeffs(P_do,s);
P_deC=coeffs(P_de,s);
S=solve(P_deC==P_doC,[wn1o,wn2o,z1o,z2o]);% Distintas soluciones para los polos propuestos
L=double(Form_Ack_obs(A,C,u_o))       % Constantes para Observabilidad 
disp('Comprobacion (comp=eig(A-L*C) ):')
comp=eig(A-L*C)     % Comprobacion 
%% Ejercicio 8
disp('EJERCICIO 8')
syms t0 tf q0 dq0 qf dqf a0 a1 a2 a3
% Polinomio cubico
disp('Polinomio cubico')
M=[1 t0 t0^2 t0^3;0 1 2*t0 3*t0^2;1 tf tf^2 tf^3;0 1 2*tf 3*tf^2];
a=[a0;a1;a2;a3];
q=[q0;dq0;qf;dqf];
a=simplify(inv(M)*q)
t0=0;   tf=2;   q0=0;   dq0=0;  dqf=1;
a=subs(a)
syms t0 tf q0 dq0 d2q0 qf dqf d2qf a0 a1 a2 a3 a4 a5
disp('Polinomio quinto grado')
M=[1 t0 t0^2 t0^3 t0^4 t0^5;...
    0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;...
    0 0 2 6*t0 12*t0^2 20*t0^3;...
    1 tf tf^2 tf^3 tf^4 tf^5;...
    0 1 2*tf 3*tf^2 4*tf^3 5*tf^4;...
    0 0 2 6*tf 12*tf^2 20*tf^3];
a=[a0;a1;a2;a3;a4;a5];
q=[q0;dq0;d2q0;qf;dqf;d2qf];
a=simplify(inv(M)*q)
t0=0;   tf=2;   q0=0;   dq0=0;  dqf=1;  d2q0=0;  d2qf=d2qf;
a=subs(a)
%% Ejercicio 11
disp('EJERCICIO 11')
disp('Primer polinomio')
syms t0 tf q0 dq0 qf dqf a0 a1 a2 a3
M=[1 t0 t0^2 t0^3;0 1 2*t0 3*t0^2;1 tf tf^2 tf^3;0 1 2*tf 3*tf^2];
a=[a0;a1;a2;a3];
q=[q0;dq0;qf;dqf];
a=simplify(inv(M)*q);
t0=0;   tf=1;   q0=0;   dq0=0;  qf=5;  dqf=40; 
a=subs(a)
disp('Segundo polinomio')
syms t0 tf q0 dq0 qf dqf a0 a1 a2 a3
M2=[1 t0 t0^2 t0^3;0 1 2*t0 3*t0^2;1 tf tf^2 tf^3;0 1 2*tf 3*tf^2];
a=[a0;a1;a2;a3];
q=[q0;dq0;qf;dqf];
a=simplify(inv(M2)*q);
t0=1;   tf=2;   q0=5;   dq0=40;  qf=10;  dqf=10; 
a=subs(a)