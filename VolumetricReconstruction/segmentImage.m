function L=segmentImage(im,epsilon)

eeF=strel('disk',round(max(size(im))/300)); 
eeG=strel('disk',5);

imS=imopen(imclose(im,eeF),eeF);
%imS=rgb2hsv(imS);
%imS=imS(:,:,1:2);

id=imdilate(imS,eeG);
ie=imerode(imS,eeG);
gr=double(squeeze(max(abs(id-ie),[],3)));
if (nargin < 2)
    epsilon=max(max(max(max(imS))))*graythresh(gr)/7;
end
plat=(gr<=epsilon);
dist=bwdist(~plat);
seeds=imregionalmax(dist);
gr(seeds>0)=-inf;
gr(gr<epsilon)=0;
L=watershed(gr);

%for l=1:max(max(max(L)))
%    props(l).size=sum(sum(sum(L==l)));
%    props(l).mean=[mean(mean(mean(d1(L==l)))) mean(mean(mean(d2(L==l))))];
%end

