clc
clear all
%% datos
load datos1.mat
t=p(:,1);
q1=q(:,2)*100;
q2=q(:,3)*100;
q3=q(:,4)*100;
p1=p(:,2)*100;
p2=p(:,3)*100;
p3=p(:,4)*100;
n=size(t);
%% Animacion
myworld=vrworld('vrml1.wrl');
open(myworld)
fig=view(myworld,'-internal');
set(myworld,'RecordInterval',[0,t(end)]);
set(fig,'Record2DFileName','Animacion.avi');
set(fig,'Record2D','on');
set(fig,'Record2DCompressQuality',100);
set(myworld,'RecordMode','scheduled');
val=0;
for i=2:n
    myworld.Auto.translation=[p1(i),0,q1(i)];
    myworld.Auto2.translation=[p2(i),0,q2(i)];
    myworld.Auto3.translation=[p3(i),0,q3(i)];
        m=(p1(i)-p1(i-1))/(q1(i)-q1(i-1));
        angulo=atan(m);
        myworld.Auto.rotation=[0,1,0,angulo-1.5708];
        myworld.Auto2.rotation=[0,1,0,angulo-1.5708];
        myworld.Auto3.rotation=[0,1,0,angulo-1.5708];    
    set(myworld,'Time',t(i));
    vrdrawnow;
end
close(myworld);
delete(myworld);



% myworld=vrworld('vrml1.wrl');
% open(myworld)
% mynodes=get(myworld,'Nodes')
% fig=view(world,'-internal');
% 
% close(myworld)


