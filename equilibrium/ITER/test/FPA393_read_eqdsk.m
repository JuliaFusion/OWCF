function [rgrid,zgrid,psirzN,Br,Bz,Bphi,rmaxis,zmaxis,rbbbs,zbbbs,rlim,zlim,B0,Ip,NH,NW] = FPA393_read_eqdsk(eqdsk_file)
% Reads eqdsk formatted equilibrium files and returns profiles required for raytracing
% jeras 26/04/2016: Modified from C:\cts\common\applic\Relray\Matlib\read_eqdsk0.m + read_eqdsk0c.m :
% Removed bugs & outdated syntax + updated for use with ITER equilibria.
%
% Input:
%  eqdsk_file : filename of eqdsk equilibrium file
% Output:
%  rgrid  : radial coordinates of grid points      (uniformly spaced)
%  zgrid  : vertical coordinates of grid points    (uniformly spaced)
%  psirzN : Normalized flux coordinates at grid points.
%  Br     : Radial magnetic field at grid points.
%  Bz     : Vertical magnetic field at grid points.
%  Bphi   : Toroidal magnetic field at grid points.
%    NOTE convention on array indexing!
%    psirzN(ir,iz) = psiN( rgrid(ir), zgrid(iz) )
%  rmaxis : Radial coordinate of magnetic axis
%  zmaxis : Vertical coordinate of magnetic axis
%  rbbbs  : Radial coordinate of plasma boundary
%  zbbbs  : Vertical coordinate of plasma boundary
%  rlim   : Radial coordinate of limiter
%  zlim   : Vertical coordinate of limiter
%  B0     : Central vacuum magnetic field
%  Ip     : Plasma current
%
fid = fopen(eqdsk_file,'r');
ccase = fscanf(fid,'%c',[8,6]); ccase=ccase';
idum = fscanf(fid,'%4i',[1,1]);
NW = fscanf(fid,'%4i',[1,1]);
NH = fscanf(fid,'%4i',[1,1]);

xdim   = fscanf(fid,'%16e',[1,1]);
zdim   = fscanf(fid,'%16e',[1,1]);
rzero  = fscanf(fid,'%16e',[1,1]);
rgrid1 = fscanf(fid,'%16e',[1,1]);
zmid   = fscanf(fid,'%16e',[1,1]);

rmaxis  = fscanf(fid,'%16e',[1,1]);
zmaxis  = fscanf(fid,'%16e',[1,1]);
ssimag  = fscanf(fid,'%16e',[1,1]);
ssibry  = fscanf(fid,'%16e',[1,1]);
bcentr  = fscanf(fid,'%16e',[1,1]);

cpasma  = fscanf(fid,'%16e',[1,1]); Ip = cpasma;
ssimag  = fscanf(fid,'%16e',[1,1]);
xdum    = fscanf(fid,'%16e',[1,1]);
rmaxis  = fscanf(fid,'%16e',[1,1]);
xdum    = fscanf(fid,'%16e',[1,1]);

zmaxis  = fscanf(fid,'%16e',[1,1]);
xdum    = fscanf(fid,'%16e',[1,1]);
ssibry  = fscanf(fid,'%16e',[1,1]);
xdum    = fscanf(fid,'%16e',[1,1]);
xdum    = fscanf(fid,'%16e',[1,1]);

f       = fscanf(fid,'%16e',[NW,1]);
pres    = fscanf(fid,'%16e',[NW,1]);
ffprim  = fscanf(fid,'%16e',[NW,1]);
pprime  = fscanf(fid,'%16e',[NW,1]);
psirz   = fscanf(fid,'%16e',[NW,NH]);
qpsi    = fscanf(fid,'%16e',[NW,1]);
nbbbs  = fscanf(fid,'%5i',[1,1]);
limitr = fscanf(fid,'%5i',[1,1]);
buf    = fscanf(fid,'%16e',[2,nbbbs]);
rbbbs = buf(1,:)';
zbbbs = buf(2,:)';
buf    = fscanf(fid,'%16e',[2,limitr]);
rlim  = buf(1,:)';
zlim  = buf(2,:)';

fclose(fid);

% Check for absurd values (>1000) in psirz and remove if necessary:
dummy = find(abs(psirz > 100));
if ~isempty(dummy)
    psirz(dummy) = 0;
end
  
NWg = [0:(NW-1)];
rgrid = rgrid1 + [0:(NW-1)]'*xdim/(NW-1);
%zgrid = zmid + ([0:(NW-1)]-(NW-1)/2)*zdim/(NW-1); % Jeras: This is a bug.
zgrid = zmid + ([0:(NH-1)]-(NH-1)/2)*zdim/(NH-1);

% jeras: Imported from read_eqdsk0c:
psil = mean(interp2(rgrid,zgrid,psirz',rbbbs,zbbbs));
psic = interp2(rgrid,zgrid,psirz',rmaxis,zmaxis);
if (psic-psil)/Ip < 0
    disp('WARNING => flux function and current direction are inconsistent' )
    disp('Changing flux function sign')
    psil = -psil;
    psic = -psic;
    psirz = -psirz;
end



dR = 0.01;
IW = 2:(NW-1);
IH = 2:(NH-1);

for nh = 1:NH
    dpsidR(IW,nh) = ...
        ( interp1(rgrid,psirz(:,nh),rgrid(IW)+dR/2,'*spline')...
        -interp1(rgrid,psirz(:,nh),rgrid(IW)-dR/2,'*spline'))/dR;
    dpsidR( 1,nh) = (interp1(rgrid,psirz(:,nh),rgrid( 1)+dR/2,'*spline')-psirz( 1,nh))/(dR/2);
    dpsidR(NW,nh) = (psirz(NW,nh)-interp1(rgrid,psirz(:,nh),rgrid(NW)-dR/2,'*spline'))/(dR/2);
end

for nw = 1:NW
    dpsidz(nw,IH) = ...
        ( interp1(zgrid,psirz(nw,:),zgrid(IH)+dR/2,'*spline')...
        -interp1(zgrid,psirz(nw,:),zgrid(IH)-dR/2,'*spline'))/dR;
    dpsidz(nw, 1) = (interp1(zgrid,psirz(nw,:),zgrid( 1)+dR/2,'*spline')-psirz(nw, 1))/(dR/2);
    dpsidz(nw,NH) = (psirz(nw,NH)-interp1(zgrid,psirz(nw,:),zgrid(NH)-dR/2,'*spline'))/(dR/2);
end

R = rgrid*ones(1,NH);   
Br = -dpsidz./(R*2*pi);                     % Use flux functions per unit angle here: normalize by 2*pi !
Bz =  dpsidR./(R*2*pi);
 
psirzN  = (psic-psirz)/(psic-psil);         % Normalize psi w.r.t value at LCFS. Adapted from read_eqdsk0c.m

fext = [f;f(NW)*ones(2*NW+1,1)];
psi  = [0:3*NW]/(NW-1);
Bphi = ones(NW,NH)*nan;
for nh = 1:NH
    Bphi(:,nh) = interp1(psi,fext,psirzN(:,nh),'*spline')./R(:,nh);
end
B0 =f(NW)/rzero;

return
