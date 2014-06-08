function vol=reconstructFromImages(clbs,data,imInd,volSize,threshold,nP,maxIt)
%clbs(i).clb = calibration
%data(i).L =
%data(i).usedL = labels of the objects to use
% nP is the number of points to use on the reconstruction process (around
% 2e5)
% maxIt is the maximum number of iterations (2-5)

global DEBUG_PLOT;
%nP=150000;
div=round(nthroot(nP,3));
minV=volSize(1:2:5)';
maxV=volSize(2:2:end)';
st=(maxV-minV)/div;
points=morePoints(minV,maxV,st);

%maxIt=5;

for iteration=1:maxIt
    if (iteration>1)
        if (size(points,2)==0)
            points=[points newPoints(nP,minV,maxV)]; %#ok<AGROW>
        else
            sP=std(points,0,2);
            minS=min(maxV-minV)/50;
            sP(sP<minS)=minS;
            coe=((maxIt-iteration)/maxIt)+1;
            minVn=(min(points,[],2)-coe*sP);
            maxVn=(max(points,[],2)+coe*sP);
            st=(maxVn-minVn)/div;
            points=[points morePoints(minVn,maxVn,st,div)]; %#ok<AGROW>
        end
    end
    
    used=ones(length(data),size(points,2),1);
    for v=1:length(data)
        proj(v).p=round(project_points(clbs(v).clb,points)); %#ok<AGROW>
        vSize(v).s=size(data(v).im); %#ok<AGROW>
    end
    
    for p=1:length(points)
        v=1;
        vals=zeros(length(data),1);

        while ((v<=length(data)))%&&(all(used(:,p))==1))
            if (    (proj(v).p(1,p)>vSize(v).s(2)) ||...
                    (proj(v).p(2,p)>vSize(v).s(1)) ||...
                    (proj(v).p(1,p)<=0)            || ...
                    (proj(v).p(2,p)<=0))
                used(v,p)=0;
            else
                vals(v)=data(v).im(proj(v).p(2,p),proj(v).p(1,p),imInd);
            end
            v=v+1;
        end
        %if all(used(p)==1)
            %if (any(vals<0.2)) %TODO REMOVE threshold
                used(vals<threshold,p)=0;
%                used(p)=0;
            %end
        %end
    end

    try
       DEBUG_PLOT; %#ok<VUNUS>
    catch %#ok<CTCH>
       DEBUG_PLOT=0;
    end
    
    if (DEBUG_PLOT==1)
        figure(1),plot3(points(1,:),points(2,:),points(3,:),'.');
        
        for ii=1:length(data)
            figure,imshow(data(ii).im(:,:,1))
            hold on
            plot(proj(ii).p(1,:),proj(ii).p(2,:),'o')
            plot(proj(ii).p(1,all(used)==1),proj(ii).p(2,all(used)==1),'.r')
            %axis tight
            
            figure(1),hold on, grid on, axis tight equal
            plot3(points(1,used(ii,:)==1),points(2,used(ii,:)==1),points(3,used(ii,:)==1),'o','Color',rand(1,3));
        end
    end
    points=points(:,all(used)==1);
end
vol=points;

    function points=morePoints(minV,maxV,st,dist)
        [a b c]=meshgrid(minV(1):st:maxV(1),minV(2):st:maxV(2),minV(3):st:maxV(3));
        a=reshape(a,[1 numel(a)]);
        b=reshape(b,[1 numel(b)]);
        c=reshape(c,[1 numel(c)]);
        points=[a;b;c];
        if (nargin>3)
            points=points+dist*(rand(size(points))-0.5);
        end
    end
    function points=newPoints(nPoints,minV,maxV)
        
        points=rand(3,nPoints);
        minV=repmat(minV,[1 size(points,2)]);
        maxV=repmat(maxV,[1 size(points,2)]);
        points=points.*(maxV-minV)+minV;
    end
end
