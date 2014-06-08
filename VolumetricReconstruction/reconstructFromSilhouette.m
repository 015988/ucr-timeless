function vol=reconstructFromSilhouette(clbs,data,volSize)
%clbs(i).clb = calibration
%data(i).L = segmented labels
%data(i).usedL = labels of the objects to use


nP=20000;

div=round(nthroot(nP,3));

minV=volSize(1:2:5)';
maxV=volSize(2:2:end)';

st=(maxV-minV)/div;

points=[];
for x=minV(1):st:maxV(1)
    for y=minV(2):st:maxV(2)
        for z=minV(3):st:maxV(3)
            points=[points [x;y;z]];
        end
    end
end


%points=newPoints(nP,minV,maxV);


for iteration=1:40
    
    if (iteration>1)
        if (size(points,2)>(nP/2))
            points=points(:,randperm(size(points,2),round(nP/2)));
        end
        
        if (size(points,2)==0)
            points=[points newPoints(nP,minV,maxV)]; %#ok<AGROW>
        else
            mP=mean(points,2);
            sP=std(points,0,2);
            sP(sP<0.1)=0.1;
            points=[points newPoints(nP-size(points,2),max(mP-10*sP,minV),min(mP+10*sP,maxV))]; %#ok<AGROW>
            points=points(:,randperm(size(points,2),nP));
        end
    end
    
    used=ones(size(points,2),1);
    
    for v=1:length(data)
        proj(v).p=round(project_points(clbs(v).clb,points)); %#ok<AGROW>
        vSize(v).s=size(data(v).L); %#ok<AGROW>
    end
    
    for p=1:length(points)
        v=1;
        while ((v<=length(data))&&(used(p)==1))
            if (    (proj(v).p(1,p)>vSize(v).s(1)) ||...
                    (proj(v).p(2,p)>vSize(v).s(2)) ||...
                    (proj(v).p(1,p)<=0)            || ...
                    (proj(v).p(2,p)<=0))
                used(p)=0;
            else
                if (~(any(data(v).L(proj(v).p(1,p),proj(v).p(2,p),1)==data(v).usedL)))
                    used(p)=0;
                end
            end
            v=v+1;
        end
    end
    
    figure(1),hold off,plot3(points(1,:),points(2,:),points(3,:),'+')
    points=points(:,used==1);
    figure(1),hold on,plot3(points(1,:),points(2,:),points(3,:),'or')
    drawnow;

   
end


for v=1:length(data)
    figure,imshow(label2rgb(data(v).L(:,:,1))),hold on,
    plot(proj(v).p(2,:),proj(v).p(1,:),'xr')
    proj(v).p=proj(v).p(:,used==1); %#ok<AGROW>
    figure,imshow(label2rgb(data(v).L(:,:,1))),hold on,
    plot(proj(v).p(2,:),proj(v).p(1,:),'xr')
end

figure,plot3(points(1,:),points(2,:),points(3,:),'+'), hold on, grid on, axis tight equal
vol=points;
    function points=newPoints(nPoints,minV,maxV)
        
        points=rand(3,nPoints);
        minV=repmat(minV,[1 size(points,2)]);
        maxV=repmat(maxV,[1 size(points,2)]);
        points=points.*(maxV-minV)+minV;
    end
end
