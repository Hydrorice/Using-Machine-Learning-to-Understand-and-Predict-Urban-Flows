clear all;close all;clc

%% simulation runtime check by plotting (u,z) and (u'w',z) 
ut_1h = ncread('453_1h.nc','uxyt');
zt_1h = ncread('453_1h.nc','zt');
upwp_1h = ncread('453_1h.nc','upwpxyt');
uw_1h = ncread('453_1h.nc','uwxyt');
ut_1h = mean(ut_1h,2);
zt_1h = mean(zt_1h,2);
upwp_1h = mean(upwp_1h,2);
uw_1h = mean(uw_1h,2);

ut_2h = ncread('453_2h.nc','uxyt');
zt_2h = ncread('453_2h.nc','zt');
upwp_2h = ncread('453_2h.nc','upwpxyt');
uw_2h = ncread('453_2h.nc','uwxyt');
ut_2h = mean(ut_2h,2);
zt_2h = mean(zt_2h,2);
upwp_2h = mean(upwp_2h,2);
uw_2h = mean(uw_2h,2);

ut_4h = ncread('453_4h.nc','uxyt');
zt_4h = ncread('453_4h.nc','zt');
upwp_4h = ncread('453_4h.nc','upwpxyt');
uw_4h = ncread('453_4h.nc','uwxyt');
ut_4h = mean(ut_4h,2);
zt_4h = mean(zt_4h,2);
upwp_4h = mean(upwp_4h,2);
uw_4h = mean(uw_4h,2);

ut_6h = ncread('xytdump.453.nc','uxyt');
zt_6h = ncread('xytdump.453.nc','zt');
upwp_6h = ncread('xytdump.453.nc','upwpxyt');
uw_6h = ncread('xytdump.453.nc','uwxyt');
ut_6h = mean(ut_6h,2);
zt_6h = mean(zt_6h,2);
upwp_6h = mean(upwp_6h,2);
uw_6h = mean(uw_6h,2);

figure

subplot(1,2,1)
plot(ut_1h,zt_1h,ut_2h,zt_2h,ut_4h,zt_4h,ut_6h,zt_6h,'linewidth',2)
title('(a) log-law check','Interpreter','latex')
xlabel("$u (m/s)$",'Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
legend('1 hr','2 hrs','4 hrs','6 hrs','Location','northwest')
ylim([0 480])
set(gca,'fontsize',25);

subplot(1,2,2)
plot(-(upwp_1h+uw_1h),zt_1h,-(upwp_2h+uw_2h),zt_2h,-(upwp_4h+uw_4h),zt_4h,-(upwp_6h+uw_6h),zt_6h,'linewidth',2)
title('(b) vertical Reynolds stresses check','Interpreter','latex')
xlabel("$-u'w' (m^2/s^2)$",'Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
legend('1 hr','2 hrs','4 hrs','6 hrs')
ylim([0 480])
set(gca,'fontsize',25);

%% domain length check by plotting (u,z) and (u'w',z) 
ut_453 = ncread('xytdump.453.nc','uxyt');
zt_453 = ncread('xytdump.453.nc','zt');
upwp_453 = ncread('xytdump.453.nc','upwpxyt');
uw_453 = ncread('xytdump.453.nc','uwxyt');
ut_453 = mean(ut_453,2);
zt_453 = mean(zt_453,2);
upwp_453 = mean(upwp_453,2);
uw_453 = mean(uw_453,2);

ut_500 = ncread('xytdump.500.nc','uxyt');
zt_500 = ncread('xytdump.500.nc','zt');
upwp_500 = ncread('xytdump.500.nc','upwpxyt');
uw_500 = ncread('xytdump.500.nc','uwxyt');
ut_500 = mean(ut_500,2);
zt_500 = mean(zt_500,2);
upwp_500 = mean(upwp_500,2);
uw_500 = mean(uw_500,2);

ut_025 = ncread('xytdump.025.nc','uxyt');
zt_025 = ncread('xytdump.025.nc','zt');
upwp_025 = ncread('xytdump.025.nc','upwpxyt');
uw_025 = ncread('xytdump.025.nc','uwxyt');
ut_025 = mean(ut_025,2);
zt_025 = mean(zt_025,2);
upwp_025 = mean(upwp_025,2);
uw_025 = mean(uw_025,2);

figure

subplot(1,2,1)
plot(ut_025,zt_025,ut_500,zt_500,ut_453,zt_453,'linewidth',2)
title('(a) log-law check','Interpreter','latex')
xlabel("$u (m/s)$",'Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
legend('025','500','453','Location','northwest')
ylim([0 480])
set(gca,'fontsize',25);

subplot(1,2,2)
plot(-(upwp_025+uw_025),zt_025,-(upwp_500+uw_500),zt_500,-(upwp_453+uw_453),zt_453,'linewidth',2)
title('(b) vertical Reynolds stresses check','Interpreter','latex')
xlabel("$-u'w' (m^2/s^2)$",'Interpreter','latex')
ylabel('height (m)','Interpreter','latex')
legend('025','500','453')
ylim([0 480])
set(gca,'fontsize',25);

%% mean vertical velocity check
vt_025 = ncread('xytdump.025.nc','vxyt');
vt_025 = mean(vt_025,'all')

vt_500 = ncread('xytdump.500.nc','vxyt');
vt_500 = mean(vt_500,'all')

vt_453 = ncread('xytdump.453.nc','vxyt');
vt_453 = mean(vt_453,'all')