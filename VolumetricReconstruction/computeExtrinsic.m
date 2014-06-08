function [cO cT err]=computeExtrinsic(A,kk,curI,curW,opt)

f_init(1)=A(1,1);
f_init(2)=A(2,2);
c_init(1)=A(1,3);
c_init(2)=A(2,3);

X_kk=[curW ; zeros(1,size(curW,2))];
xn = normalize_pixel(curI,f_init,c_init);
Np = size(xn,2);

% Check for planarity of the structure:
X_mean = mean(X_kk,2);
R_transform = eye(3) ;
T_transform = -(R_transform)*X_mean;
X_new = R_transform*X_kk + T_transform*ones(1,Np);
% Compute the planar homography:
H = compute_homography(xn,X_new(1:2,:));
% De-embed the motion parameters from the homography:
sc = mean([norm(H(:,1));norm(H(:,2))]);
H = H/sc;

u1 = H(:,1);
u1 = u1 / norm(u1);
u2 = H(:,2) - dot(u1,H(:,2)) * u1;
u2 = u2 / norm(u2);
u3 = cross(u1,u2);
RRR = [u1 u2 u3];
Rckk = rodrigues(rodrigues(RRR));
Tckk = H(:,3);

%If Xc = Rckk * X_new + Tckk, then Xc = Rckk * R_transform * X_kk + Tckk + T_transform
cT = Tckk + Rckk* T_transform;
Rckk = Rckk * R_transform;
cO = rodrigues(Rckk);
if (opt==1)
    global DEBUG_PLOT; %#ok<TLEV>
    DEBUG_PLOT=0;
    [nX,err]=fminsearch(@fcn_opt_OT,[cO; cT],optimset('Display','off','TolFun',10,'MaxFunEvals',10E2,'MaxIter',10E2));
    if (err<10E3)
        [nX,err]=fminsearch(@fcn_opt_OT,nX);%'TolFun',1,'MaxFunEvals',10E3,'MaxIter',10E3));
    end
    cO=nX(1:3);
    cT=nX(4:6);
end
    function vErr=fcn_opt_OT(xIn)
        vErr=frameError(curI, curW,...
            [A(1,1) A(2,2) A(1,3) A(2,3) kk(1) kk(2)], xIn);
    end

end