function xOut=project_points(calib,xIn)
%function xOut=projection(calib,xIn)
%xOut: (u,v) -> projection of xIn by calib
%xIn: 3D xIn world points
%calib: calibration (load_calibration)

if (min(size(xIn))>=3)
    if (size(xIn,1)>size(xIn,2))
        xIn=xIn';
    end
end
xIn(4,:)=1;

if (isstruct(calib))
    RT=calib.RT;
    KK=calib.KK;
    k=calib.dir.k;
    p=calib.dir.p;
    s=calib.dir.s;
else
    [KK RT k p s]=convert_vet_calib(calib);
end

point=RT*xIn;
point=point./repmat(point(3,:),[3 1]);

u=point(1,:);
v=point(2,:);
%distortion

u2v2=u.^2+v.^2;

du=k(1).*u.*u2v2+k(2).*u.*(u2v2).^2+s(1).*u2v2+(p(1).*(3.*u.^2+v.^2)+2.*p(2).*u.*v);
dv=k(1).*v.*u2v2+k(2).*v.*(u2v2).^2+s(2).*u2v2+(p(2).*(3.*v.^2+u.^2)+2.*p(1).*u.*v);
point(1,:)=u+du; %( 4334.34 4038802 )
point(2,:)=v+dv;

point=KK*point;
point=point./repmat(point(3,:),[3 1]);

xOut=point(1:2,:);
