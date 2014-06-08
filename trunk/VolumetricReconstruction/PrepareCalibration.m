function clb=PrepareCalibration(I,W,A,kk,statCam)
%function [A OT usedFrames]=PrepareCalibration(W,I,imHeight,imWidth)
%Performs *extrinsic* camera calibration using the world points W, image points I.
% statCam: if the camera is static (1) or not (0)
% if KK is provided, extrinsic parameters are computed
% if KK is *not* provied, A is computed, but only A.
% size(W)=(2*nFrames,nPoints), size(I)=(2*nFrames,nPoints)
% Based on a heavily modified version of
% init_intrinsic_param.m/compute_extrinsic.m from J. Bouguet
% http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html

global DEBUG_PLOT;
DEBUG_PLOT=0;

realFrames=find((sum(I(1:2:end,:)==-1,2)+sum(W(1:2:end,:)==-1,2))==0);
if (~statCam)
    for p=1:size(I,2)
        x=I(1:2:end,p);
        xl=interp1(find(x>-1),x(x>-1),1:length(x),'spline');
        I(1:2:end,p)=xl;
        
        y=I(2:2:end,p);
        yl=interp1(find(y>-1),y(y>-1),1:length(y),'spline');
        I(2:2:end,p)=yl;
    end
    
    for p=1:size(W,2)
        x=W(1:2:end,p);
        xl=interp1(find(x>-1),x(x>-1),1:length(x));
        W(1:2:end,p)=xl;
        
        y=W(2:2:end,p);
        yl=interp1(find(y>-1),y(y>-1),1:length(y));
        W(2:2:end,p)=yl;
    end
end
usedFrames=find((sum(I(1:2:end,:)==-1,2)+sum(W(1:2:end,:)==-1,2))==0);
nFrames=size(I,1)/2;

%% estimating extrinsic parameters

T=zeros(nFrames,3);
O=zeros(nFrames,4);
for f=1:nFrames
    if (isempty(find(usedFrames==f,1)))
        T(f,:)=nan(1,3);
        O(f,:)=nan(1,4);
    else

        curI=squeeze(I(((2*(f-1))+1):(2*f),:));
        curW=squeeze(W(((2*(f-1))+1):(2*f),:));
        [tempO T(f,:) err]=computeExtrinsic(A,kk,curI,curW,1);
        if (err<5) %orig 100
            O(f,:)=dcm2quat(rodrigues(tempO));
        else
            usedFrames=setdiff(usedFrames,f);
            T(f,:)=nan(1,3);
            O(f,:)=nan(1,4);
        end
    end
end

if (~statCam)
    for p=1:3
        tmp=T(:,p);
        if (any(isnan(tmp)))
            idx=find(~isnan(tmp));
            tmp=interp1(idx,tmp(~isnan(tmp)),1:nFrames,'cubic')';
            T(:,p)=tmp;
        end
    end
    
    for p=1:4
        tmp=O(:,p);
        if (any(isnan(tmp)))
            idx=find(~isnan(tmp));
            tmp=interp1(idx,tmp(~isnan(tmp)),1:nFrames,'cubic')';
            O(:,p)=tmp;
        end
    end
else
    %there is probably a better way of doing this bit
    for p=1:3
        tmp=T(:,p);
        if (any(isnan(tmp)))

            tmp(isnan(tmp))=nanmean(tmp);
            T(:,p)=tmp;
        end
    end
    
    for p=1:4
        tmp=O(:,p);
        if (any(isnan(tmp)))
            tmp(isnan(tmp))=nanmean(tmp);
            O(:,p)=tmp;
        end
    end
    
end

for f=1:nFrames
    clb(f).clb.RT=[quat2dcm(O(f,:)) T(f,:)'];
    clb(f).clb.KK=A;
    clb(f).clb.dir.k=kk;
    clb(f).clb.dir.p=[0 0];
    clb(f).clb.dir.s=[0 0];
end

    function Aerr=fcnOptAOT(x)
        cIn=x(1:6);
        cOT=reshape(x(7:end),sOT);
        Aerr=0;
        for cf=usedFrames'
            curI=squeeze(I(((2*(cf-1))+1):(2*cf),:));
            curW=squeeze(W(((2*(cf-1))+1):(2*cf),:));
            Aerr=Aerr+frameError(curI,curW,cIn,cOT(cf,:));
        end
    end
end

