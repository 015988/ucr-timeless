function [trajRet spTraj]=trajDefBall(allPoints,MaxFunEval,step)
TT=sort(unique(allPoints(:,5)));
cP=1;
z=zeros(1,length(TT));
x=z;
y=z;
for t=TT'
    %curPoints=allPoints(allPoints(:,5)==t,1:3);
    temp=median(allPoints(allPoints(:,5)==t,1:3));
    %figure,plot3(curPoints(:,1),curPoints(:,2),curPoints(:,3),'.')
    %hold on, grid on, axis tight equal
    x(cP)=temp(1);
    y(cP)=temp(2);
    z(cP)=temp(3);
    cP=cP+1;
end
%step=10; %the borders can be somewhat twitchy
trajOrig=[[x(2) x(3:step:(end-2)) x(end-1)]; [y(2) y(3:step:(end-2)) y(end-1)]; [z(2) z(3:step:(end-2)) z(end-1)]; [TT(2) TT(3:step:(end-2))' TT(end-1)]];
sT=size(trajOrig);
x0=reshape(trajOrig,[1 numel(trajOrig)]);
if (nargin==1)
    t=fminsearch(@fcn_opt,x0,optimset('TolFun',10,'Display','final','TypicalX',x0));%,'MaxFunEval',1000));
else
    t=fminsearch(@fcn_opt,x0,optimset('TolFun',10,'Display','final','TypicalX',x0,'MaxFunEval',MaxFunEval));
end
traj=reshape(t,sT);
%traj=sortrows(traj',4)';
spTraj=csapi(traj(4,:),traj(1:3,:));
trajRet=fnval(spTraj,TT);%spTraj.breaks(1):splDiv:spTraj.breaks(end))';
trajRet=[trajRet; TT'];

uPos=allPoints(:,4)>0;
uNeg=allPoints(:,4)<=0;
figure,fnplt(spTraj),hold on,
plot3(allPoints(uPos,1),allPoints(uPos,2),allPoints(uPos,3),'.g','MarkerSize',2);
plot3(allPoints(uNeg,1),allPoints(uNeg,2),allPoints(uNeg,3),'.r','MarkerSize',1);
hold on, grid on, axis tight equal

    function err=fcn_opt(x)
        err=trajFitBall(reshape(x,sT),25,allPoints); %diameter of the ball 50mm
    end

end