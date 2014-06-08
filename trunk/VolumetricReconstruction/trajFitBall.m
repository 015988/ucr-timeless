function err=trajFitBall(trajPoints,r,points)

splDiv=0.5;
%trajPoints=sortrows(trajPoints',4)';
spTraj=csapi(trajPoints(4,:),trajPoints(1:3,:));
%figure(9),hold off, 
%fnplt(spTraj),hold on, 
%plot3(trajPoints(1,:),trajPoints(2,:),trajPoints(3,:),'or');
%plot3(points(:,1),points(:,2),points(:,3),'.g','MarkerSize',1);
%grid on, axis tight equal

% tic,idx=rangesearch(traj,points(:,1:3),r);toc
% 
% used=zeros(length(idx),1);
% for i=1:length(idx)
%     if (~isempty(idx{i}))
%         used(i)=1;
%     end
% end

traj=fnval(spTraj,spTraj.breaks(1):splDiv:spTraj.breaks(end))';

d=pdist2(traj,points(:,1:3));
used=sum(d<=r,1)>0;

err=-sum(points(used==1,4));
%title(sprintf('err= %g',err));
%pause(0.5)
end