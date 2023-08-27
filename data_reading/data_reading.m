clear all; close all; clc;

% reading file
file = 'tdump.560.nc';

% % obtain the description of netCDF file
% ncdisp(file)

%% log-law check
% ut = ncread(file,'ut');% 'Streamwise velocity'
% ut = mean(ut,4);
% 
% zt = ncread(file,'zt');% 'Vertical displacement of cell centers'
% 
% fig=figure;
% 
%     % method 1
%     ut_m=mean(ut);
%     ut_m=mean(ut_m);
%     ut_m=squeeze(ut_m);
%     subplot(3,1,1)
%     plot(ut_m,zt)
%     
%     % method 2
%     utt=zeros(1,length(zt));
%     for i=1:length(zt)
%         utt(i)=mean2(ut(:,:,i));
%     end
%     subplot(3,1,2)
%     plot(utt,zt)
% 
%     % method 3 (need xytdump file)
% %     ut = ncread('xytdump.025.nc','uxyt');
% %     zt = ncread('xytdump.025.nc','zt');
% %     subplot(3,1,3)
% %     plot(ut,zt)
%     
% han=axes(fig,'visible','off'); 
% han.Title.Visible='on';
% han.XLabel.Visible='on';
% han.YLabel.Visible='on';
% ylabel(han,'vertical displacement z (m)');
% xlabel(han,'streamwise velocity u (m/s)');
% title(han,'log-law');
 
%% u'w' check
% upwp = ncread('xytdump.453.nc','upwpxyt');
% upwp = mean(upwp,2);
% 
% uw = ncread('xytdump.453.nc','uwxyt');
% uw = mean(uw,2);
% 
% ut = ncread('xytdump.453.nc','uxyt');
% ut = mean(ut,2);
% 
% zt = ncread('xytdump.453.nc','zt');
% 
% figure
% subplot(1,2,1)
% plot(-(upwp+uw) ,zt,'linewidth',2)
% title('(a)')
% xlabel("$u'w' (m^2/s^2)$",'Interpreter','latex')
% ylabel('height (m)','Interpreter','latex')
% ylim([0 480])
% set(gca,'fontsize',25);
% 
% % from figure, (u_star)^2=0.56
u_star = 0.75;
% 
% 
% subplot(1,2,2)
% plot(ut/u_star,zt,ones(size(zt)),zt,'--','linewidth',2)
% title('(b)')
% xlabel("$u/u^\star$",'Interpreter','latex')
% ylabel('height (m)','Interpreter','latex')
% ylim([0 480])
% set(gca,'fontsize',25);
% 
% % from figure, delta=63
delta = 63;

%% y-z plane flow profile check for flow shifting (only for driver simulation)

% ut = ncread(file,'ut');
% ut = mean(ut,4);
% 
% xt = ncread(file,'xt');
% yt = ncread(file,'yt');
% zt = ncread(file,'zt');
% 
% utc = zeros(length(xt),length(yt),length(zt));
% for i=1:length(xt)-1 
%     utc(i,:,:) = (ut(i,:,:) + ut(i+1,:,:)) / 2; 
% end
% utc(end,:,:) = ut(end,:,:);
% 
% x_i= length(xt)/2;
% utc_b = squeeze(utc(1,:,:));
% utc_m = squeeze(utc(x_i,:,:));
% utc_l = squeeze(utc(end,:,:));
% [Y,Z] = meshgrid(zt,yt);
% 
% figure
% subplot(3,1,1)
% contour(Z,Y,utc_b)
% title('the streamwise velocity profile in the Y-Z plane at the beginning of the domain','Interpreter','latex')
% xlabel('y (m)','Interpreter','latex');ylabel('z (m)','Interpreter','latex')
% ylim([0,500]);
% set(gca,'fontsize',15);
% colorbar
% 
% subplot(3,1,2)
% contour(Z,Y,utc_m)
% title('the streamwise velocity profile in the Y-Z plane at the mid-length of the domain','Interpreter','latex')
% xlabel('y (m)','Interpreter','latex');ylabel('z (m)','Interpreter','latex')
% ylim([0,500]);
% set(gca,'fontsize',15);
% colorbar
% 
% subplot(3,1,3)
% contour(Z,Y,utc_l)
% title('the streamwise velocity profile in the Y-Z plane at the end of the domain','Interpreter','latex')
% xlabel('y (m)','Interpreter','latex');ylabel('z (m)','Interpreter','latex')
% ylim([0,500]);
% set(gca,'fontsize',15);
% colorbar

%% x-y plane flow profile check for flow shifting (only for driver simulation)

% load('bluered.mat')
% 
% ut = ncread(file,'ut');
% ut = mean(ut,4);
% 
% vt = ncread(file,'vt');
% vt = mean(vt,4);
% 
% z_index = 70;
% 
% ut_xy=ut(:,:,z_index);
% vt_xy=vt(:,:,z_index);
% 
% figure
% 
% subplot(2,1,1)
% pcolor(ut_xy')
% title('the streamwise velocity in the X-Y plane','Interpreter','latex')
% xlabel('grid number in the x direction','Interpreter','latex')
% ylabel('grid number in the y direction','Interpreter','latex')
% set(gca,'fontsize',15);
% axis equal
% colorbar
% 
% subplot(2,1,2)
% pcolor(vt_xy')
% title('the spanwise velocity in the X-Y plane','Interpreter','latex')
% xlabel('grid number in the x direction','Interpreter','latex')
% ylabel('grid number in the y direction','Interpreter','latex')
% set(gca,'fontsize',15);
% axis equal
% colorbar
% 
% ut_mean = mean(ut,'all')
% vt_mean = mean(vt,'all')

%% x,y,z coordinates

xt = ncread(file,'xt');% 'West-East displacement of cell centers'
xm = ncread(file,'xm');% 'West-East displacement of cell edges'

yt = ncread(file,'yt');% 'South-North displacement of cell centers'
ym = ncread(file,'ym');% 'South-North displacement of cell edges'

zt = ncread(file,'zt');% 'Vertical displacement of cell centers'
zm = ncread(file,'zm');% 'Vertical displacement of cell edges'

time = ncread(file,'time');% 'tstatsdump, Output time for statistics'

%% velocities

ut = ncread(file,'ut');% 'Streamwise velocity'
ut = mean(ut,4);

vt = ncread(file,'vt');% 'Spanwise velocity'
vt = mean(vt,4);

wt = ncread(file,'wt');% 'Vertical velocity'
wt = mean(wt,4);

pt = ncread(file,'pt');% 'Pressure'
pt = mean(pt,4);

%% velocities at cell center 

% x-direction
utc = zeros(length(xt),length(yt),length(zt));
for i=1:length(xt)-1 
    utc(i,:,:) = (ut(i,:,:) + ut(i+1,:,:)) / 2; 
end
utc(end,:,:) = ut(end,:,:);

% y-direction
vtc = zeros(length(xt),length(yt),length(zt));
for j=1:length(yt)-1 
    vtc(:,j,:) = (vt(:,j,:) + vt(:,j+1,:)) / 2; 
end
vtc(:,end,:) = vt(:,end,:);

% z-direction
wtc = zeros(length(xt),length(yt),length(zt));
for k=1:length(zt)-1 
    wtc(:,:,k) = (wt(:,:,k) + wt(:,:,k+1)) / 2;  
end
wtc(:,:,end) = wt(:,:,end);

%% Reynolds Stresses (Turbulent momentum fluxes)

uu = ncread(file,'upuptc');% 'u variance'
uu = mean(uu,4);

vv = ncread(file,'vpvptc');% 'v variance'
vv = mean(vv,4);

ww = ncread(file,'wpwptc');% 'w variance'
ww = mean(ww,4);

uv = ncread(file,'upvpt');% 'uv co-variance'
uv = mean(uv,4);

uw = ncread(file,'upwpt');% 'uw co-variance'
uw = mean(uw,4);

vw = ncread(file,'vpwpt');% 'vw co-variance'
vw = mean(vw,4);

%% Finite Difference method

% rename grid coordinate
x = xt;
y = yt;
z = zt;

% grid size
dx = x(2)-x(1);
dy = y(2)-y(1);
dz = z(2)-z(1);

% grid numbers
nx = length(x);
ny = length(y);
nz = length(z);

% extend the datasets (ghost points)
% with respect to the x-direction
ut_x_left = utc(1,:,:)-(utc(2,:,:)-utc(1,:,:));
ut_x_right= utc(end,:,:)+utc(end,:,:)-utc(end-1,:,:);

pt_x_left = pt(1,:,:)-(pt(2,:,:)-pt(1,:,:));
pt_x_right= pt(end,:,:)+pt(end,:,:)-pt(end-1,:,:);

uu_x_left = uu(1,:,:)-(uu(2,:,:)-uu(1,:,:));
uu_x_right= uu(end,:,:)+uu(end,:,:)-uu(end-1,:,:);

% with respect to the y-direction
ut_y_forward = utc(:,1,:)-(utc(:,2,:)-utc(:,1,:));
ut_y_backward= utc(:,end,:)+utc(:,end,:)-utc(:,end-1,:);

uv_y_forward = uv(:,1,:)-(uv(:,2,:)-uv(:,1,:));
uv_y_backward= uv(:,end,:)+uv(:,end,:)-uv(:,end-1,:);

% with respect to the z-direction
ut_z_bottom = utc(:,:,1)-(utc(:,:,2)-utc(:,:,1));
ut_z_top= utc(:,:,end)+utc(:,:,end)-utc(:,:,end-1);

uw_z_bottom = uw(:,:,1)-(uw(:,:,2)-uw(:,:,1));
uw_z_top= uw(:,:,end)+uw(:,:,end)-uw(:,:,end-1);

%% Backward Finite Difference

% approximation of derivative d/dx
dudx = zeros(nx,ny,nz);
dudx(1,:,:) = (utc(1,:,:)-ut_x_left)/dx;

dpdx = zeros(nx,ny,nz);
dpdx(1,:,:) = (pt(1,:,:)-pt_x_left)/dx;

duudx = zeros(nx,ny,nz);
duudx(1,:,:) = (uu(1,:,:)-uu_x_left)/dx;

for i = 2:nx
    dudx(i,:,:) = (utc(i,:,:)-utc(i-1,:,:))/dx;
    dpdx(i,:,:) = (pt(i,:,:)-pt(i-1,:,:))/dx;
    duudx(i,:,:) = (uu(i,:,:)-uu(i-1,:,:))/dx;
end

% approximation of derivative d/dy
dudy = zeros(nx,ny,nz);
dudy(:,1,:) = (utc(:,1,:)-ut_y_forward)/dy;

duvdy = zeros(nx,ny,nz);
duvdy(:,1,:) = (uv(:,1,:)-uv_y_forward)/dy;

for j = 2:ny
    dudy(:,j,:) = (utc(:,j,:)-utc(:,j-1,:))/dy;
    duvdy(:,j,:) = (uv(:,j,:)-uv(:,j-1,:))/dy;
end

% approximation of derivative d/dz
dudz = zeros(nx,ny,nz);
dudz(:,:,1) = (utc(:,:,1)-ut_z_bottom)/dy;

duwdz = zeros(nx,ny,nz);
duwdz(:,:,1) = (uw(:,:,1)-uw_z_bottom)/dy;

for k = 2:nz
    dudz(:,:,k) = (utc(:,:,k)-utc(:,:,k-1))/dy;
    duwdz(:,:,k) = (uw(:,:,k)-uw(:,:,k-1))/dy;
end

%% laplacian & 2nd order derivative by central FD

% extend datasets
ut_ex = zeros(nx+2,ny,nz);
for i = 1:nx
    ut_ex(i+1,:,:) = utc(i,:,:);
end
ut_ex(1,:,:) = ut_x_left;
ut_ex(end,:,:) = ut_x_right;

ut_ey = zeros(nx,ny+2,nz);
for j = 1:ny
    ut_ey(:,j+1,:) = utc(:,j,:);
end
ut_ey(:,1,:) = ut_y_forward;
ut_ey(:,end,:)=ut_y_backward;

ut_ez = zeros(nx,ny,nz+2);
for k = 1:nz
    ut_ez(:,:,k+1) = utc(:,:,k);
end
ut_ez(:,:,1) = ut_z_bottom;
ut_ez(:,:,end)=ut_z_top;

% 2nd order derivative by central FD
d2udx2 = zeros(nx,ny,nz);
for i = 2:nx+1
    d2udx2(i-1,:,:)= (ut_ex(i+1,:,:)-2*ut_ex(i,:,:)+ut_ex(i-1,:,:))/dx^2;
end

d2udy2 = zeros(nx,ny,nz);
for j = 2:ny+1
    d2udy2(:,j-1,:)= (ut_ey(:,j+1,:)-2*ut_ey(:,j,:)+ut_ey(:,j-1,:))/dy^2;
end

d2udz2 = zeros(nx,ny,nz);
for k = 2:nz+1
    d2udz2(:,:,k-1)= (ut_ez(:,:,k+1)-2*ut_ez(:,:,k)+ut_ez(:,:,k-1))/dz^2;
end

% laplacian
lap_u = d2udx2 + d2udy2 + d2udz2;

%% side view at 322.5 m with y-index equal to 65

y_i = 65;% y-index

% dimension reduction
utc_sv = squeeze(utc(:,y_i,:));
vtc_sv = squeeze(vtc(:,y_i,:));
wtc_sv = squeeze(wtc(:,y_i,:));
dpdx_sv = squeeze(dpdx(:,y_i,:));
dudx_sv = squeeze(dudx(:,y_i,:));
dudy_sv = squeeze(dudy(:,y_i,:));
dudz_sv = squeeze(dudz(:,y_i,:));
duudx_sv = squeeze(duudx(:,y_i,:));
duvdy_sv = squeeze(duvdy(:,y_i,:));
duwdz_sv = squeeze(duwdz(:,y_i,:)); 
lap_u_sv = squeeze(lap_u(:,y_i,:));

% dominant balance model
uux_sv = utc_sv.*dudx_sv;
vuy_sv = vtc_sv.*dudy_sv;
uuxvuy_sv = uux_sv + vuy_sv;% first term

wuz_sv = wtc_sv.*dudz_sv;% second term

uvpyuupx_sv = duvdy_sv + duudx_sv;% fifth term

%% plotting side view

% rip out the data of buildings 

% ordinary buildings
for i = 0:17
    N = 21+i*32;
    uuxvuy_sv(N:N+7,1:33) = 0;
    wuz_sv(N:N+7,1:33) = 0;
    dpdx_sv(N:N+7,1:33) = 0;
    uvpyuupx_sv(N:N+7,1:33) = 0;
    duwdz_sv(N:N+7,1:33) = 0;
    lap_u_sv(N:N+7,1:33) = 0;
end

% the high-rising building
uuxvuy_sv(117:124,1:97) = 0;
wuz_sv(117:124,1:97) = 0;
dpdx_sv(117:124,1:97) = 0;
uvpyuupx_sv(117:124,1:97) = 0;
duwdz_sv(117:124,1:97) = 0;
lap_u_sv(117:124,1:97) = 0;

% load distinguished colour map
load('bluered.mat')

% 2-D grid coordinates generation
[X,Z] = meshgrid(x,z);

% kinematic viscosity, due to the use of LES, nu is near zero and could be 
% neglected in RANS.
nu = 1.5e-5;

figure('position',[150 175 1210 480]);
sgtitle('side view at y = 322.5 m','Interpreter','latex','fontsize',15)


subplot(3,2,1)

% original figure
% pcolor(X',Z',uuxvuy_sv)
% title('$\bar{u}\bar{u}_x+\bar{v}\bar{u}_y$','Interpreter','latex')

% normalized figure
pcolor(X',Z',uuxvuy_sv/(u_star^2/delta))
title('$(\bar{u}\bar{u}_x+\bar{v}\bar{u}_y)/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
ylim([0 400])
%caxis([-0.1 0.1])
caxis([-10 10])

subplot(3,2,2)

% original figure
% pcolor(X',Z',wuz_sv)
% title('$\bar{w}\bar{u}_z$','Interpreter','latex')

% normalized figure
pcolor(X',Z',wuz_sv/(u_star^2/delta))
title('$\bar{w}\bar{u}_z/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
ylim([0 400])
% caxis([-0.1 0.1])
caxis([-10 10])

subplot(3,2,3)

% original figure
% pcolor(X',Z',dpdx_sv)
% title('$\rho^{-1}\bar{p}_x$','Interpreter','latex')

% normalized figure
pcolor(X',Z',dpdx_sv/(u_star^2/delta))
title('$\bar{p}_x/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
ylim([0 400])
% caxis([-0.1 0.1])
caxis([-10 10])


% subplot(3,2,4)
% pcolor(X',Z',nu*lap_u_sv)
% title('$\nu\nabla^2\bar{u}$','Interpreter','latex')
% xlabel('length (m)','Interpreter','latex')
% ylabel('height (m)','Interpreter','latex')
% shading flat
% colorbar()
% colormap(cmap)
% axis equal tight
% %xlim([0 1500])
% caxis([-0.1 0.1])


subplot(3,2,4)

% original figure
% pcolor(X',Z',uvpyuupx_sv)
% title('$\overline{{u^\prime}^2}_x+\overline{u^\prime v^\prime}_y$','Interpreter','latex')

% normalized figure
pcolor(X',Z',uvpyuupx_sv/(u_star^2/delta))
title('$(\overline{{u^\prime}^2}_x+\overline{u^\prime v^\prime}_y)/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
ylim([0 400])
% caxis([-0.1 0.1])
caxis([-10 10])

subplot(3,2,5)

% original figure
% pcolor(X',Z',duwdz_sv)
% title('$\overline{u^\prime w^\prime}_z$','Interpreter','latex')

% normalized figure
pcolor(X',Z',duwdz_sv/(u_star^2/delta))
title('$\overline{u^\prime w^\prime}_z/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
ylim([0 400])
% caxis([-0.1 0.1])
caxis([-10 10])

% plot the sum of RANS terms
rans_sv = uuxvuy_sv+wuz_sv+dpdx_sv+uvpyuupx_sv+duwdz_sv;
rans_sv_mean = mean(rans_sv,'all','omitnan');

subplot(3,2,6)

% original figure
pcolor(X',Z',rans_sv)
title('the sum of RANS terms','Interpreter','latex')

% normalized figure
pcolor(X',Z',rans_sv/(u_star^2/delta))
title('the sum of RANS terms/$\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
ylim([0 400])
% caxis([-0.1 0.1])
caxis([-10 10])

%% plan view at 148.75 m with z-index equal to 60

% z_i = 20;% z-index, plan view at 48.75 m
% z_i = 40;% z-index, plan view at 98.75 m
z_i = 60;% z-index

% dimension reduction
utc_pv = squeeze(utc(:,:,z_i));
vtc_pv = squeeze(vtc(:,:,z_i));
wtc_pv = squeeze(wtc(:,:,z_i));
dpdx_pv = squeeze(dpdx(:,:,z_i));
dudx_pv = squeeze(dudx(:,:,z_i));
dudy_pv = squeeze(dudy(:,:,z_i));
dudz_pv = squeeze(dudz(:,:,z_i));
duudx_pv = squeeze(duudx(:,:,z_i));
duvdy_pv = squeeze(duvdy(:,:,z_i));
duwdz_pv = squeeze(duwdz(:,:,z_i)); 
lap_u_pv = squeeze(lap_u(:,:,z_i));

% dominant balance model
uuxvuy_pv = utc_pv.*dudx_pv + vtc_pv.*dudy_pv;% first term

wuz_pv = wtc_pv.*dudz_pv;% second term

uvpyuupx_pv = duvdy_pv + duudx_pv;% fifth term

%% plotting plane view

% rip out the data of buildings 

% the high-rising building
uuxvuy_pv(117:124,61:68) = 0;
wuz_pv(117:124,61:68) = 0;
dpdx_pv(117:124,61:68) = 0;
uvpyuupx_pv(117:124,61:68) = 0;
duwdz_pv(117:124,61:68) = 0;
lap_u_pv(117:124,61:68) = 0;

% 2-D grid coordinates generation
[X,Y] = meshgrid(x,y);

figure('position',[150 175 1210 480]);
sgtitle('plan view at z = 148.75 m','Interpreter','latex','fontsize',15)


subplot(3,2,1)

% original figure
% pcolor(X',Y',uuxvuy_pv)
% title('$\bar{u}\bar{u}_x+\bar{v}\bar{u}_y$','Interpreter','latex')

% normalized figure
pcolor(X',Y',uuxvuy_pv/(u_star^2/delta))
title('$(\bar{u}\bar{u}_x+\bar{v}\bar{u}_y)/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('width (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
% caxis([-0.1 0.1])
caxis([-10 10])

subplot(3,2,2)

% original figure
% pcolor(X',Y',wuz_pv)
% title('$\bar{w}\bar{u}_z$','Interpreter','latex')

% normalized figure
pcolor(X',Y',wuz_pv/(u_star^2/delta))
title('$\bar{w}\bar{u}_z/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('width (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
% caxis([-0.1 0.1])
caxis([-10 10])

subplot(3,2,3)

% original figure
% pcolor(X',Y',dpdx_pv)
% title('$\rho^{-1}\bar{p}_x$','Interpreter','latex')

% normalized figure
pcolor(X',Y',dpdx_pv/(u_star^2/delta))
title('$\bar{p}_x/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('width (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
% caxis([-0.1 0.1])
caxis([-10 10])


% subplot(3,2,4)
% pcolor(X',Y',nu*lap_u_pv)
% title('$\nu\nabla^2\bar{u}$','Interpreter','latex')
% xlabel('length (m)','Interpreter','latex')
% ylabel('width (m)','Interpreter','latex')
% shading flat
% colorbar()
% colormap(cmap)
% axis equal tight
% %xlim([0 1500])
% caxis([-0.1 0.1])


subplot(3,2,4)

% original figure
% pcolor(X',Y',uvpyuupx_pv)
% title('$\overline{{u^\prime}^2}_x+\overline{u^\prime v^\prime}_y$','Interpreter','latex')

% normalized figure
pcolor(X',Y',uvpyuupx_pv/(u_star^2/delta))
title('$(\overline{{u^\prime}^2}_x+\overline{u^\prime v^\prime}_y)/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('width (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
% caxis([-0.1 0.1])
caxis([-10 10])


subplot(3,2,5)

% original figure
% pcolor(X',Y',duwdz_pv)
% title('$\overline{u^\prime w^\prime}_z$','Interpreter','latex')

% normalized figure
pcolor(X',Y',duwdz_pv/(u_star^2/delta))
title('$\overline{u^\prime w^\prime}_z/\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('width (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
% caxis([-0.1 0.1])
caxis([-10 10])

% plot the sum of RANS terms
rans_pv = uuxvuy_pv+wuz_pv+dpdx_pv+uvpyuupx_pv+duwdz_pv;
rans_pv_mean = mean(rans_pv,'all','omitnan');


subplot(3,2,6)

% original figure
% pcolor(X',Y',rans_pv)
% title('the sum of RANS terms','Interpreter','latex')

% normalized figure
pcolor(X',Y',rans_pv/(u_star^2/delta))
title('the sum of RANS terms/$\frac{u^\star{^2}}{\delta}$','Interpreter','latex')

xlabel('length (m)','Interpreter','latex')
ylabel('width (m)','Interpreter','latex')
set(gca,'fontsize',15);
shading flat
colorbar()
colormap(cmap)
axis equal tight
xlim([0 1500])
% caxis([-0.1 0.1])
caxis([-10 10])
