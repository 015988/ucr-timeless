function [L, seq]=segmentSequence(vidName,frameRange,epsilon)

v=VideoReader(vidName);

if (nargin==1)
    frameRange=[1 v.NumberOfFrames];
end

seq=read(v,frameRange);

 for i=1:size(seq,4)
     tmp=rgb2hsv(squeeze(seq(:,:,:,i)));
     nD(:,:,i)=tmp(:,:,2);
 end

 seq=nD;
 clear nD
 
ee=strel('disk',9,6);
seq=imopen(imclose(seq,ee),ee); 


id=imdilate(seq,ones(3,3,3));
ie=imerode(seq,ones(3,3,3));


gr=double(squeeze(max(abs(id-ie),[],4)));
%figure,imshow(gr(:,:,1),[0 200])
%gr=imclose(gr,ones(7,7,2));%figure,imshow(gr(:,:,1),[0 200])
%gr=imopen(gr,ones(3,3,1));%figure,imshow(gr(:,:,1),[0 200])


if (nargin <= 2)
    epsilon=max(max(max(max(seq))))*graythresh(gr)/6;
end

plat=gr<=epsilon;
dist=bwdist(~plat);
seeds=imregionalmax(dist);
gr(seeds>0)=-inf;
gr(gr<epsilon)=0;
L=watershed(gr);



