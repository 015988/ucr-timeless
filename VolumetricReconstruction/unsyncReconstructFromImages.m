function points=unsyncReconstructFromImages(clbs,data,volSize,threshold,nP,minPoints,maxIt,refVid,outRange)
%clbs(i).clb = calibration
%refVidId is a number, representing which video to use as reference
% nP is the number of points to use on the reconstruction process
% maxIt is the maximum number of iterations (2-5)
% range in which points with < thr still are preserved (outside range)

global DEBUG_PLOT;

try
    DEBUG_PLOT; %#ok<VUNUS>
catch %#ok<CTCH>
    DEBUG_PLOT=0;
end
if (isempty(DEBUG_PLOT))
    DEBUG_PLOT=0;
end

div=round(nthroot(5*nP,3));
minV=volSize(1:2:5)';
vMinLast=minV;
maxV=volSize(2:2:end)';
vMaxLast=maxV;
st=(maxV-minV)/div;
points=morePoints(minV,maxV,st);
notFoundMult=1;

iteration=1;
while (iteration<=maxIt)
    if ((iteration>1)||(notFoundMult>1))
        if (size(points,2)==0)
            points=[points newPoints(round(notFoundMult*5*nP),vMinLast,vMaxLast)]; %#ok<AGROW>
        else
            sP=std(points(1:3,:),0,2);
            minS=min(maxV-minV)/50;
            sP(sP<minS)=minS;
            coe=((maxIt-iteration)/maxIt)+1;
            minVn=max(vMinLast,(min(points(1:3,:),[],2)-coe*sP));
            maxVn=min(vMaxLast,(max(points(1:3,:),[],2)+coe*sP));
            st=(maxVn-minVn)/div;
            points=[points morePoints(minVn,maxVn,st,div)]; %#ok<AGROW>
        end
    end
    
    vmax=-ones(length(data),size(points,2),1);
    
    for v=1:length(data)
        vSize(v).s=size(data(v).im); %#ok<AGROW>
        for f=1:size(data(v).im,3)
            
            tmpProj=round(project_points(clbs(v).clb(f).clb,points(1:3,:)));
            used=ones(size(points,2),1);
            used((tmpProj(1,:)>vSize(v).s(2)))=0;
            used((tmpProj(2,:)>vSize(v).s(1)))=0;
            used((tmpProj(1,:)<=0))=0;
            used((tmpProj(2,:)<=0))=0;
            
            tmpProj=tmpProj(:,used==1);
            points=points(:,used==1);
            vmax=vmax(:,used==1);
            
            tInd=sub2ind(size(data(v).im),tmpProj(2,:),tmpProj(1,:),repmat(f,[1 size(tmpProj,2)]));
            vmax(v,:)=max(vmax(v,:),data(v).im(tInd));
            if (v==refVid)
                points(5,(data(v).im(tInd)>0)&(vmax(v,:)==data(v).im(tInd)))=f;
            end
        end
    end
    
    [used points(4,:)]=selectPoints(vmax,threshold);
    
    vMinLast=max(vMinLast,min(points(1:3,:),[],2));
    vMaxLast=min(vMaxLast,max(points(1:3,:),[],2));
    
    outPoints=points(:,used==0);
    if (size(outPoints,2)>nP) %if we dont keep outPoints small, pdist2 will eat the whole memory
        outPoints=outPoints(:,randperm(size(outPoints,2),nP));
    end
    points=points(:,used==1);
    
    
    if (size(points,2)>=minPoints)
        iteration=iteration+1;
        if (iteration>maxIt)
            outPoints(4,:)=-1;
            outPoints=outPoints(:,min(pdist2(outPoints(1:3,:)',points(1:3,:)'),[],2)<outRange);
            %figure,plot3(points(1,:),points(2,:),points(5,:),'.b'); hold on, grid on
            %plot3(outPoints(1,:),outPoints(2,:),outPoints(5,:),'.r')
            if (~isempty(points))
                points=[points outPoints]; %#ok<AGROW>
            end
        end
    else
        if (notFoundMult>16)
            return
        else
            notFoundMult=1.2*notFoundMult;
        end
    end
end

    function points=morePoints(minV,maxV,st,dist)
        [a b c]=meshgrid(minV(1):st:maxV(1),minV(2):st:maxV(2),minV(3):st:maxV(3));
        a=reshape(a,[1 numel(a)]);
        b=reshape(b,[1 numel(b)]);
        c=reshape(c,[1 numel(c)]);
        points=[a;b;c];
        if (nargin>3)
            points=points+dist*(rand(size(points))-0.5);
        end
        points(4,:)=0;
        points(5,:)=0;
    end
    function points=newPoints(nPoints,minV,maxV)
        
        points=rand(3,nPoints);
        minV=repmat(minV,[1 size(points,2)]);
        maxV=repmat(maxV,[1 size(points,2)]);
        points=points.*(maxV-minV)+minV;
        points(4,:)=0;
        points(5,:)=0;
    end
    function [used vals]= selectPoints(vmax,threshold)
        %vmax=sort(vmax);
        %used=mean(vmax(1:3,:))>threshold;
        %used=median(vmax)>=threshold;
        used=mean(vmax)>=threshold;
        vals=mean(vmax);
        %used=vmax>=threshold;
        %used=(sum(used)==size(used,1))|(sum(used)>=(size(used,1)-1) & (min(vmax,[],1)<=(threshold/3)));
    end
end
