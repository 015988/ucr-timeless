function exLine

W=[0:0.2:1; 0:0.2:1; 0:0.2:1];

figure(9),plot3(W(1,:),W(2,:),W(3,:),'.'),hold on, grid on, axis tight equal


A1=[100 0 400;0 100 800;0 0 1];
T1=-[5 -2 2]'; 
R1=eye(3); 
%u1=generateProjection(A1,R1,T1,W(:,1:end));
%figure,plot(u1(1,:),u1(2,:),'x-');

A2=[100 0 400;0 100 800;0 0 1];
R2=eye(3);
T2=-[-2 5 2]'; 
%u2=generateProjection(A2,R2,T2,W(:,1:end));
%figure,plot(u2(1,:),u2(2,:),'o-');

intersPlanes(W,-T1,'r');
intersPlanes(W,-T2,'b');
axis tight equal