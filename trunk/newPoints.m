function points=newPoints(nPoints,minV,maxV)

points=rand(3,nPoints);
minV=repmat(minV,[1 size(points,2)]);
maxV=repmat(maxV,[1 size(points,2)]);
points=points.*(maxV-minV)+minV;