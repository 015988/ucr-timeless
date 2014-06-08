function [A kk]=IntrinsicCalibration(I,W,imHeight,imWidth,useDistortion)
%function [A OT usedFrames]=IntrinsicCalibration(W,I,imHeight,imWidth)
%Performs camera calibration using the world points W, image points I.
% imHeight and imWidth == size of the original
%image (used to initialize principal point).
% useDistortion = 0/1
% size(W)=(2*nFrames,nPoints), size(I)=(2*nFrames,nPoints)
% Based on a heavily modified version of
% init_intrinsic_param.m/compute_extrinsic.m from J. Bouguet
% http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html

global DEBUG_PLOT;
DEBUG_PLOT=0;

c_init = [imWidth;imHeight]/2 - 0.5; % initialize at the center of the image
% matrix that subtract the principal point:
Sub_cc = [1 0 -c_init(1);0 1 -c_init(2);0 0 1];

% Compute explicitely the focal length using all the (mutually orthogonal) vanishing points
% note: The vanishing points are hidden in the planar collineations H_kk

usedFrames=find((sum(I(1:2:end,:)==-1,2)+sum(W(1:2:end,:)==-1,2))==0);
nValids=length(usedFrames);
nFrames=size(I,1)/2;

V=zeros(2*nValids,2);
b=zeros(2*nValids,1);
usedLines=1;
for f=1:nFrames
    if (isempty(find(usedFrames==f,1)))
        cal(f).H=[]; %#ok<*AGROW>
    else
        curI=squeeze(I(((2*(f-1))+1):(2*f),:));
        curW=squeeze(W(((2*(f-1))+1):(2*f),:));
        cal(f).H=compute_homography(curI,curW);
        [A_kk b_kk]=extractVanishingPoints(cal(f).H,Sub_cc);
        V((2*(usedLines-1)+1):(2*usedLines),:)=A_kk;
        b((2*(usedLines-1)+1):(2*usedLines),:)=b_kk;
        usedLines=usedLines+1;
    end
end

f_init = sqrt(abs(1./(((V'*V)\(V'))*b)));
% Initial estimate for intrinsic parameters
A = [f_init(1) 0 c_init(1);0 f_init(2) c_init(2); 0 0 1];

%small bundle adjustment to get A right
allOT=zeros(nValids,6);
corresp=zeros(nValids,1);
frErrors=zeros(nValids,1);
cF=1;
for f=usedFrames'
    curI=squeeze(I(((2*(f-1))+1):(2*f),:));
    curW=squeeze(W(((2*(f-1))+1):(2*f),:));
    [o T]=computeExtrinsic(A,[0 0],curI,curW,0);
    frErrors(cF)=frameError(curI,curW,[f_init(1) f_init(2) c_init(1) c_init(2) 0 0],[o; T]);
    allOT(cF,:)=[o;T];
    corresp(cF)=f;
    cF=cF+1;
end

% how many frames (in decreasing order of error)
k_clusters=min(20,nFrames);
IDX=kmeans(allOT,k_clusters,'emptyaction','drop');
selected_idx=[];
for i=1:k_clusters
    tErr=frErrors;
    tErr(IDX~=i)=NaN;
    [~,cIdx]=max(tErr);
    if (~isempty(cIdx))
        selected_idx=[selected_idx;corresp(cIdx)];
    end
end
usedFrames=sort(selected_idx);
nValids=length(usedFrames);
OT=zeros(nValids,6);

%get only the interesting frames (aka bigger error)
cF=1;
for f=usedFrames'
    curI=squeeze(I(((2*(f-1))+1):(2*f),:));
    curW=squeeze(W(((2*(f-1))+1):(2*f),:));
    [o T]=computeExtrinsic(A,[0 0], curI,curW,1);
    OT(cF,:)=[o;T];
    cF=cF+1;
end

sOT=size(OT);

if (useDistortion==1)
    x0=[0 0];
    %sprintf('Calibration error (pre): %g',fcnkk(x0))
    kk=fminsearch(@fcnkk,x0);
    %sprintf('Calibration error (pos): %g',fcnkk(kk))

    x0=[ A(1,1) A(2,2) A(1,3) A(2,3) kk(1) kk(2) reshape(OT,[1 numel(OT)])];
    [xS cErr]=fminsearch(@fcnOptAkk,x0,optimset('TypicalX',x0+0.1,'Display','off'));
    kk=xS(5:6);
else
    x0=[ A(1,1) A(2,2) A(1,3) A(2,3) reshape(OT,[1 numel(OT)])];
    [xS cErr]=fminsearch(@fcnOptA,x0,optimset('TypicalX',x0+0.1,'Display','off'));
    kk=[0 0];

end

A=[xS(1) 0    xS(3);...
    0   xS(2) xS(4);...
    0    0     1];

    function Aerr=fcnkk(x)
        cIn=[A(1,1) A(2,2) A(1,3) A(2,3) x(1) x(2)];
        cOT=OT;
        Aerr=0;
        nf=1;
        for cf=usedFrames'
            curI=squeeze(I(((2*(cf-1))+1):(2*cf),:));
            curW=squeeze(W(((2*(cf-1))+1):(2*cf),:));
            Aerr=Aerr+frameError(curI,curW,cIn,cOT(nf,:));
            nf=nf+1;
        end
    end


    function Aerr=fcnOptAkk(x)
        cIn=x(1:6);
        cOT=reshape(x(7:end),sOT);
        Aerr=0;
        nf=1;
        for cf=usedFrames'
            curI=squeeze(I(((2*(cf-1))+1):(2*cf),:));
            curW=squeeze(W(((2*(cf-1))+1):(2*cf),:));
            Aerr=Aerr+frameError(curI,curW,cIn,cOT(nf,:));
            nf=nf+1;
        end

    end
    function Aerr=fcnOptA(x)
        cIn=x(1:4);
        cIn(5)=0;
        cIn(6)=0;
        cOT=reshape(x(5:end),sOT);
        Aerr=0;
        nf=1;
        for cf=usedFrames'
            curI=squeeze(I(((2*(cf-1))+1):(2*cf),:));
            curW=squeeze(W(((2*(cf-1))+1):(2*cf),:));
            Aerr=Aerr+frameError(curI,curW,cIn,cOT(nf,:));
            nf=nf+1;
        end

    end


end

