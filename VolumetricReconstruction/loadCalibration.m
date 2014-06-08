function calib=loadCalibration(name)
%function calib=load_calibration(name)
%name=filename of calibration. If '', uses pick('clb')
%Returns a structure with:
%P: projection matrix (DLT)
%KK: 3x3 intrinsic matrix
%c: camera center
%p : decentering
%s : thin prism
%k: radial distortion
%
%Further details:
%J. Wang, F. Shia, J. Zhanga, and Y. Liu. A new calibration model of camera lens distortion. Pattern
%Recognition, 41(2):607–615, February 2008.
%(we do not use the proposed method, but the explained in the intro)

if (~exist('radial','var'))
    radial=1;
end

%empty name == pick
if (isempty(name))
    name=pick('clb');
end

% try to open the file
f=fopen(name,'r');

%loading saved results
K=fscanf(f,'%e',9);
RT=fscanf(f,'%e',12);

calib.KK=unstack_matrix(K,3,3);
calib.RT=unstack_matrix(RT,3,4);

calib.R=calib.RT(1:3,1:3);
calib.T=calib.RT(:,4);
calib.P=calib.KK*calib.RT;
calib.dir.k=fscanf(f,'%e %e\n',2);
calib.dir.p=fscanf(f,'%e %e\n',2);
try
    calib.dir.s=fscanf(f,'%e %e\n',2);
    
    calib.inv.k=fscanf(f,'%e %e\n',2);
    calib.inv.p=fscanf(f,'%e %e\n',2);
    calib.inv.s=fscanf(f,'%e %e\n',2);
catch %#ok<CTCH>
    %if it had a failure reading the further values, it's from an older
    %version of the toolbox.
    calib.inv.k=calib.dir.p;
    calib.dir.p=[0;0];
    calib.dir.s=[0;0];
    calib.inv.p=[0;0];
    calib.inv.s=[0;0];
end
fclose(f);
%additional calculations
calib.KKi=pinv(calib.KK);
calib.Ri=pinv(calib.R);
calib.c=-calib.Ri*calib.KKi*calib.P(:,4);
disp(sprintf('Camera %s loaded',name));
end