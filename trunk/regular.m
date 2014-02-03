function regular

minV=[ 1  1  1]'-1;
maxV=[11 11 11]'+1;


rng(9)
W=generate3d(6,300,minV'+1,maxV'-1);
figure(9),plot3(W(1,:),W(2,:),W(3,:),'.'),hold on, grid on, axis tight equal
figure(10),plot(sqrt(sum((W(:,1:(size(W,2)-1))-W(:,2:end)).^2))*100,'+-'),hold on, grid on %speed

A1=[100 0 400;0 100 800;0 0 1];
T1=-[30 5 20]'; 
R1=eye(3); 
u1=generateProjection(A1,R1,T1,W(:,1:3:end));

A2=[100 0 400;0 100 800;0 0 1];
R2=eye(3);
T2=-[5 30 20]'; 
u2=generateProjection(A2,R2,T2,W(:,2:3:end));

A3=[200 0 200;0 200 200;0 0 1];
R3=eye(3); 
T3=-[-20 -20 20]';
u3=generateProjection(A3,R3,T3,W(:,3:3:end));


for i=1:size(u1,2)
    cu1=u1(:,i);
    cu2=u2(:,i);
    cu3=u3(:,i);
    if (i==1)
        xinit=[0 0 0];
    else
        xinit=res(i-1,:);
    end
    t=tic;
    [res(i,:) fval,exitflag,output]=fminsearch(@fcn_opt,xinit,optimset('MaxIter',5E3,'TolX',1E-16,'TolFun',1E-16,'MaxFunEvals',1E12));%
    toc(t);
end

W1=W(:,1:3:end);
W2=W(:,2:3:end);
W3=W(:,3:3:end);

er1=sqrt(sum((W1-res').^2));
er2=sqrt(sum((W2-res').^2));
er3=sqrt(sum((W3-res').^2));
figure,plot(er1,'x-'),hold on, grid on
print('-dpng','-r600','reg_W1_err.png')
figure,plot(er2,'x-'),hold on, grid on
print('-dpng','-r600','reg_W2_err.png')
figure,plot(er3,'x-'),hold on, grid on
print('-dpng','-r600','reg_W3_err.png')

save regular3Cam
dir

    function err=fcn_opt(x)
        x1=generateProjection(A1,R1,T1,x');
        x2=generateProjection(A2,R2,T2,x');
        x3=generateProjection(A3,R3,T3,x');
        err=sqrt(sum((cu1-x1).^2))+...
            sqrt(sum((cu2-x2).^2))+...
            sqrt(sum((cu3-x3).^2));
    end
        

end