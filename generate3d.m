function data=generate3d(nControlPoints,nPoints,minW,maxW)
base=rand(3,nControlPoints);
minW=repmat(minW,[nControlPoints 1])';
maxW=repmat(maxW,[nControlPoints 1])';
base=base.*(maxW-minW)+minW;
%figure,plot3(base(1,:),base(2,:),base(3,:),'+-')
%hold on, grid on
data=csapi(0:(nControlPoints-1),base,[0:(nControlPoints-1)/(nPoints):(nControlPoints-1) nControlPoints-1]);
data=data(:,1:nPoints);
%plot3(data(1,:),data(2,:),data(3,:),'+-r')
end