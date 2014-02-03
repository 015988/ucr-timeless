function [nUsed timeNear]=pointsToCurveDist(points,spCurve,epsilon,mask,getTime)
%timeNear=pointsToCurveDist(points,spCurve,epsilon)
% nUsed = points used
% timeNear is a cell, same size of nUsed, with empty value if the point is
% more distant than epsilon of the curve, or 't' such that spCurve(t) is
% the closest point.
% if isempty(mask), process all points, otherwise only the ones in mask
% if getTime~=1, it does not properly refines the involved time.

if (isempty(mask))
    nUsed=[];
    timeNear=[];
    return;
else
    points=points(:,mask);
end

splDiv=0.001;
dist=fnval(spCurve,spCurve.breaks(1):splDiv:spCurve.breaks(end));
timeNear=rangesearch(dist',points',epsilon);

nUsed=zeros(size(points,2),1);

parfor i=1:length(mask)
    if (~isempty(timeNear{i}))
        if (getTime==1)
            [tTemp fval]=fminbnd(@(t)sum(sqrt(sum((fnval(spCurve,t)-points(:,i)).^2))),...
                max(spCurve.breaks(1),mean(timeNear{i}*splDiv)-1),...
                min(spCurve.breaks(end),mean(timeNear{i}*splDiv)+1),...
                optimset('TolX',1e-2,'UseParallel','never','MaxIter',50));
            if (fval<epsilon)
                timeNear{i}=tTemp;
                nUsed(i)=1;
            else
                timeNear{i}=[];
            end
        else
            nUsed(i)=1;
        end
    end
end

nUsed=find(nUsed)';
if (getTime==1)
    newNear=[];
    for i=1:length(nUsed)
        newNear{i}=timeNear{nUsed(i)}; %#ok<AGROW>
    end
    timeNear=newNear;
end

for i=1:length(nUsed)
    nUsed(i)=mask(nUsed(i));
end




