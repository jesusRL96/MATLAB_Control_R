clc
close all
clear all
% ejercicio 7:
% +++++++++++++++ DATOS A INGRESAR +++++++++++++++
% a,alpha,d,theta
M={'l1','-pi/2','0','pi+q1';...  %Datos obtenidos por medio de Denavit Hartenberg
    'l2','0','0','q2';...
    'l3','0','0','q3'};
syms g
G=[0;g;0];  %direccion de la gravedad
J=[1,1,1];  %Tipo de junta: 1-->Revoluta; 0-->Prismatica

% ++++++++++++++++++ RESULTADOS ++++++++++++++++++
%% Matrices A
disp('                * Matrices A')
[f,c]=size(M);
As={};
for i=1:f
    v = genvarname(strcat('A',num2str(i))); 
    eval([v ' = CalculoA(M(i,:))']);
end
% Multiplicaciones de matrices A
disp('                * Multiplicaciones de matrices A')
for i=1:f
    v = genvarname(strcat('A',num2str(i),'s')); 
end
si='1';
for i=2:f
    vact = genvarname(strcat('A',num2str(i)));
    vant = genvarname(strcat('A',si));
    si=strcat(si,num2str(i));
    v = genvarname(strcat('A',si));
    eval([v,'=simplify(',vant,'*',vact,')'])
end
%% Z_n, O_n y Ocm_n
disp('                * Valores de Z_n, O_n y Ocm_n')
Z_s=[0;0;1];    %canonicos
O_s=[0;0;0];
si='';
for i=1:f-1
    %Z_n
    si=strcat(si,num2str(i));
    v = genvarname(strcat('A',si));
    eval([strcat('Z_s=[Z_s,',v,'(1:3,3)','];')]);
    %O_n
    eval([strcat('O_s=[O_s,',v,'(1:3,4)','];')]);
end
Ocm_s=[A1(1:3,4)];
Ocm_s=subs(Ocm_s,'l1','lc1');
si='1';

for i=2:f
    si=strcat(si,num2str(i));
    v = genvarname(strcat('A',si));
    l = strcat('l',num2str(i));
    lc = strcat('lc',num2str(i));
    
    %Ocm_n 
    eval([strcat('O = ',v,'(1:3,4)',';')]);
    eval([strcat('O = ','subs(','O,''',l,''',''',lc,''');')])
    Ocm_s=[Ocm_s,O];
end
Z_s
O_s
Ocm_s
%% Jacobianos
disp('                *Jacobianos:')
for i=1:f
    Jcm=sym(zeros(6,f));
    
    if J(i)==1      %Revoluta
        for j=1:i
            Jv=simplify(cross(Z_s(:,j),Ocm_s(:,i)-O_s(:,j)));
            Jw=Z_s(:,j);
            Jcm(:,j)=[Jv;Jw];
        end
        
    else            %Prismatica
        for j=1:i
            Jv=Z_s(:,j);
            Jw=zeros(3,1);
            Jcm(:,j)=[Jv;Jw];
        end
    end
    v = genvarname(strcat('Jcm',num2str(i)));
    disp(strcat('            -Jacobiano:',num2str(i)))
    eval([strcat(v,'=Jcm')]);
end
%% Matriz D(q)
disp('             *Matriz de inercia D(q)')
m_s=sym(zeros(f,1));
for i=1:f       %Matrices de Inercia Diagonales
    Ix=sym(strcat('Ix',num2str(i)));
    Iy=sym(strcat('Iy',num2str(i)));
    Iz=sym(strcat('Iz',num2str(i)));
    v = genvarname(strcat('I',num2str(i)));
    eval([strcat(v,'=[Ix,0,0;0,Iy,0;0,0,Iz];')]);
    v = genvarname(strcat('m',num2str(i)));
    mi=eval([strcat('sym(v)')]);
    m_s(i)=mi;
end
m_s;

Jvi=[];
Jwi=[];
Ri=[];
D=0;
si='';
for i=1:f       %Matriz D(q)
    si=strcat(si,num2str(i));
    Di=0;
    v = genvarname(strcat('Jcm',num2str(i)));
    Jvi=eval([strcat(v,'(1:3,:)')]);
    Jwi=eval([strcat(v,'(4:end,:)')]);
    v = genvarname(strcat('A',si));
    Ri=eval([strcat(v,'(1:3,1:3)')]);
    v = genvarname(strcat('I',num2str(i)));
    Ii=eval([v]);
    mi=m_s(i);
    Di=simplify((mi*transpose(Jvi)*Jvi) + (transpose(Jwi)*Ri*Ii*transpose(Ri)*Jwi));
    D=simplify(D + Di);
end    
D=sym(D);
D=expand(D);
D=simplify(D)
%% Matriz C(q,qp)
disp('                  *Matriz de coriolis C(q,qp)')
[f,c]=size(D);
q=[];
qp=[];
syms qi
for i=1:f
    v = genvarname(strcat('q',num2str(i)));
    qi=eval([strcat('sym(v);')]);
    q=[q,qi]; 
    v = genvarname(strcat('q',num2str(i),'p'));
    qip=eval([strcat('sym(v);')]);
    qp=[qp,qip];   
end

for i=1:f
    for j=1:f
        for k=1:f
            dkj = D(k,j);
            dki = D(k,i);
            dij = D(i,j);
            v = genvarname(strcat('C',num2str(i),num2str(j),num2str(k)));
            Cijk=((diff(dkj,q(i),1) + diff(dki,q(j),1) - diff(dij,q(k),1))/2);
            Cijk=simplify(Cijk);
            eval([strcat(v,'=Cijk;')])           
        end
    end
end
C=zeros(f);
C=sym(C);
for k=1:f
    for j=1:f
        for i=1:f
            v = genvarname(strcat('C',num2str(i),num2str(j),num2str(k)));
            Cijk = eval([v]);
            Ckj = Cijk*qp(i);
            C(k,j) = C(k,j)+ Ckj;
            C(k,j) = simplify(C(k,j));
        end
    end
end
C=sym(C);
C=expand(C);
C=simplify(C)
%% Matriz g(q)
disp('                     *vector de efectos gravitacionales g(q)')
rci_s=sym(zeros(3,f));
rci_s(2,:)=Ocm_s(2,:);
P=0;
for i=1:f
    v = genvarname(strcat('m',num2str(i)));
    mi=m_s(i);
    rci=rci_s(:,i);
    Pi=mi*transpose(G)*rci;
    P=P+Pi;
end
P=simplify(P);
phi=[];
for i=1:f
    qi=q(i);
    phi_i=diff(P,qi,1);
    phi=[phi;phi_i];
end
phi=simplify(phi)
%% fuerzas
disp('                *Vector de fuerzas T')
T=[];
for i=1:f
    if J(i)==1
        v = genvarname(strcat('T',num2str(i)));
        Ti=eval([strcat('sym(v);')]);
    else
        v = genvarname(strcat('F',num2str(i)));
        Ti=eval([strcat('sym(v);')]);
    end
        T=[T;Ti];
end
T
%% Valores numericos
disp('VALORES NUMERICOS')
l1=.3;  l2=.25; l3=.16;
lc1=.15;    lc2=.125;   lc3=.08;
Ix1=.02;    Iy1=.012;   Iz1=.01;
Ix2=.01;    Iy2=.012;   Iz2=.01;
Ix3=.08;    Iy3=.018;   Iz3=.01;
m1=9;   m2=5;   m3=3;
D=vpa(subs(D))
C=vpa(subs(C))
phi=vpa(subs(phi))
%% Propiedad de antisimetria
disp('                  *PROPIEDAD DE ANTISIMETRIA')
S=0;
for i=1:f
    Si=simplify(diff(D,q(i))*qp(i));
    S=simplify(S+Si);
end
disp('S=Dp-2*C:')
S=S-(2*C);
S=sym(S);
S=expand(S);
S=simplify(S)
disp('S+S'':')
S+transpose(S)
