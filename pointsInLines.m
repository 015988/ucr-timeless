function newPoints=pointsInLines(points)



newPoints=zeros(3,3*size(points,2));

pointCorrespBetween=randperm(size(points,2));
pointCorrespLine=randperm(size(points,2));
for i=1:size(points,2)
    newPoints(:,i)=(points(:,i)+points(:,pointCorrespBetween(i)))/2;
    newPoints(:,(  size(points,2)+i))= (points(:,i)-points(:,pointCorrespLine(i))) + points(:,i); 
    newPoints(:,(2*size(points,2)+i))= (points(:,i)-points(:,pointCorrespLine(i))) + points(:,pointCorrespLine(i));
end





end