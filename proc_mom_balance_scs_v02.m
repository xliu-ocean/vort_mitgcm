%% get the INTEGRATED momentum balance along an isopycnal surface
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Dec. 8, 2021 created

ichk = 0;

% reading grid information
xc = rdmds('grid/XC');
yc = rdmds('grid/YC');
xg = rdmds('grid/XG');
yg = rdmds('grid/YG');
rc = rdmds('grid/RC');
dxc = rdmds('grid/DXC');
dyc = rdmds('grid/DYC');
dxg = rdmds('grid/DXG');
dyg = rdmds('grid/DYG');
rAw = rdmds('grid/RAW');
raz = rdmds('grid/RAZ');
rAc = rdmds('grid/RAC');
drF = rdmds('grid/DRF');
drC = rdmds('grid/DRC');
hFacW = rdmds('grid/hFacW');
hFacS = rdmds('grid/hFacS');
% calc coriolis
om0 = 2*pi/(86400);
% f = 2*om0*sind(yc);
g = 9.81;

[nx,ny] = size(xc);
nz = size(drF,3);

f = zeros(nx,ny);
f(:,1) = 1.e-4;
for i = 2:ny
   f(:,i) = f(:,i-1) + 1.e-11*(yc(:,i)-yc(:,i-1));
end

beta = zeros(nx,ny);
beta(:,2:end) = (f(:,2:end)-f(:,1:end-1))./dyc(:,2:end);
rf(1) = 0;
for ik = 1:56
    rf(ik+1) = rf(ik)-drF(ik);
end
%% read terms
disp('Reading data from files ...')
for im = 1:1 %% for 1 years 
    disp( ['Record   ',num2str(im),' ...'])
    istr = num2str(im*5184+12374208,'%010d');
    mom_terms0 = rdmds(['diags/mom_trend.',istr]);
    dis_terms0 = rdmds(['diags/diss.',istr]);
    rho_terms0 = rdmds(['diags/state_3d.',istr]);
    eta_terms0 = rdmds(['diags/xyMom_2d.',istr]);
    if im==1
        mom_terms_m = mom_terms0;
        dis_terms_m = dis_terms0;
        rho_terms_m = rho_terms0;
        eta_terms_m = eta_terms0;
    else
        mom_terms_m = cat(5,mom_terms_m,mom_terms0);
        dis_terms_m = cat(5,dis_terms_m,dis_terms0);
        rho_terms_m = cat(5,rho_terms_m,rho_terms0);
        eta_terms_m = cat(5,eta_terms_m,eta_terms0);
    end
end

for im = 1:1
    disp( ['Record   ',num2str(im),' ...'])
    istr = num2str(im*5184+12374208,'%010d');
    xymom_terms0 = rdmds(['diags/xyMom_ave.',istr]);
    if im==1
        xymom_terms_m = xymom_terms0;
    else
        xymom_terms_m = cat(5,xymom_terms_m,xymom_terms0);
    end
end
xymom_terms = mean(xymom_terms_m,5);
mom_terms = mean(mom_terms_m,5);
eta_terms = mean(eta_terms_m,5);
dis_terms = mean(dis_terms_m,5);
rho_terms = mean(rho_terms_m,5);

% get every variable
uvel_diag = xymom_terms(:,:,:,1);
vvel_diag = xymom_terms(:,:,:,2);
wvel_diag = xymom_terms(:,:,:,3);
rhoa = squeeze(rho_terms(:,:,:,1));
theta = squeeze(rho_terms(:,:,:,2));
salt = squeeze(rho_terms(:,:,:,3));
eta = squeeze(eta_terms_m(:,:,1));

nt = size(mom_terms,5);
%% name data using the momentum terms
%'Um_Diss ' 'Vm_Diss ' 'Um_Advec' 'Vm_Advec' 'Um_Cori ' 'Vm_Cori ' 'Um_dPHdx' 
%'Vm_dPHdy' 'Um_Ext  ' 'Vm_Ext  ' 
% **** NO use ['Um_AdvZ3' 'Vm_AdvZ3' 'Um_AdvRe' 'Vm_AdvRe']
% in mom_trend 
Um_Diss = squeeze(mom_terms(:,:,:,1,:));
Vm_Diss = squeeze(mom_terms(:,:,:,2,:));
Um_Advec = squeeze(mom_terms(:,:,:,3,:)); %%(coriolis included)
Vm_Advec = squeeze(mom_terms(:,:,:,4,:)); %%(coriolis included)
Um_Cori = squeeze(mom_terms(:,:,:,5,:));
Vm_Cori = squeeze(mom_terms(:,:,:,6,:));
Um_dPH = squeeze(mom_terms(:,:,:,7,:)); %% baroclinic prs grd
Vm_dPH = squeeze(mom_terms(:,:,:,8,:)); %% baroclinic prs grd
Um_Ext = squeeze(mom_terms(:,:,:,9,:));
Vm_Ext = squeeze(mom_terms(:,:,:,10,:));

Um_Adv_a = Um_Advec - Um_Cori;  % get advection seperately
Vm_Adv_a = Vm_Advec - Vm_Cori;  % advection only

VISrI_Um0 = squeeze(dis_terms(:,:,:,4,:));
VISrI_Vm0 = squeeze(dis_terms(:,:,:,8,:));
VISCx_Um = squeeze(dis_terms(:,:,:,1,:));
VISCx_Vm = squeeze(dis_terms(:,:,:,5,:));
VISCy_Um = squeeze(dis_terms(:,:,:,2,:));
VISCy_Vm = squeeze(dis_terms(:,:,:,6,:));
VISrE_Um0 = squeeze(dis_terms(:,:,:,3,:));
VISrE_Vm0 = squeeze(dis_terms(:,:,:,7,:));
% get surface pressure grad
% for barotropic pressure
etan = squeeze(eta_terms(:,:,1,:));

Um_PHsur0 = zeros(nx,ny,nt);  Vm_PHsur0 = zeros(nx,ny,nt);  
Um_PHsur0(2:end,:,:) = -g*(etan(2:end,:,:)-etan(1:end-1,:,:))./dxc(2:end,:);
Vm_PHsur0(:,2:end,:) = -g*(etan(:,2:end,:)-etan(:,1:end-1,:))./dyc(:,2:end);
Um_PHsur = [];Vm_PHsur = [];
for ik=1:nz
    Um_PHsur(:,:,ik,:) = Um_PHsur0;
    Vm_PHsur(:,:,ik,:) = Vm_PHsur0;
end
Um_PH_tot = Um_dPH + Um_PHsur;
Vm_PH_tot = Vm_dPH + Vm_PHsur;

% get implicit vertical viscosity terms
% tendency from this vertical flux:
% Um_Impl(:,:,k) =[ VISrI_Um(:,:,k+1) - VISrI_Um(:,:,k) ]
%                 /[ rAw(:,:)*drF(k)*hFacW(:,:,k) ]
VISrI_Um = cat(3,VISrI_Um0,VISrI_Um0(:,:,1,:)*0); %% (set boundary to zeros)
VISrI_Vm = cat(3,VISrI_Vm0,VISrI_Vm0(:,:,1,:)*0); %% (set boundary to zeros)

Um_Impl = (VISrI_Um(:,:,2:end,:)-VISrI_Um(:,:,1:end-1,:))./(rAw.*drF.*hFacW);
Vm_Impl = (VISrI_Vm(:,:,2:end,:)-VISrI_Vm(:,:,1:end-1,:))./(rAw.*drF.*hFacS);
% get total of the terms
Um_tot = Um_Diss+Um_Advec+Um_dPH+Um_Ext+Um_PHsur+Um_Impl;
Vm_tot = Vm_Diss+Vm_Advec+Vm_dPH+Vm_Ext+Vm_PHsur+Vm_Impl;
if ichk
    k=40
    figure 
    plot(Um_tot(:,40,k),'k','linewidth',2)
    hold on
    plot(Um_Advec(:,40,k),'r')
    plot(Um_dPHdx(:,40,k))
    plot(Um_Ext(:,40,k))
%     plot(Um_PHsur(:,40,k))
    plot(Um_Impl(:,40,k))
    plot(Um_Diss(:,40,k),'g','linewidth',2)
end
%% get integrated terms of the momentum equations
% 0. find the depth of isopycnals

% \\\\ set the referenced diapycnal here \\\\\\\\\\
% sigref = 4.3;
sigref = 5.5;
%%%%

talpha = 2.E-4;
 sbeta = 7.4E-4;
tref = 20;
sref = 30;
rhoNil = 999.8;
rhoCon = 999.8;
dRho = rhoNil-rhoCon;
% copied from find_rho.F (MITgcm/model/src/find_rho.F)
%         rhoLoc(i,j)=rhoNil*(
%     &     sBeta*(sFld(i,j)-refSalt)
%     &   -tAlpha*(tFld(i,j)-refTemp) )
%     &        + dRho
rhoLoc = rhoNil.*(sbeta.*(salt-sref)-talpha.*(theta-tref))+dRho;
rhoLoc(rhoLoc<0) = nan;

rhoLoc_u = rhoLoc; rhoLoc_v = rhoLoc;
rhoLoc_u(2:end,:,:,:) = 0.5*(rhoLoc(1:end-1,:,:,:)+rhoLoc(2:end,:,:,:));
rhoLoc_v(:,2:end,:,:) = 0.5*(rhoLoc(:,1:end-1,:,:)+rhoLoc(:,2:end,:,:));

% get depth of isopycnals
for it = 1:nt
    disp(['Step 0: Searching thermocline, time=  ',num2str(it)])
    
    [k0,k0u,k0v] = get_layer_ref1(rhoLoc(:,:,:,it),rhoLoc_u(:,:,:,it),...
                                            rhoLoc_v(:,:,:,it),sigref);
    kref(:,:,it) = k0;                  
    kref_u(:,:,it) = k0u;                  
    kref_v(:,:,it) = k0v;
% find the position in the grid (hsec)
    [hfrac0,hfracf0,kref_fc0] = get_layer_frac(kref(:,:,:,it),rhoLoc(:,:,:,it),rc,rf,sigref);
    [hfracu0,hfracfu0,kref_fcu0] = get_layer_frac(kref(:,:,:,it),rhoLoc(:,:,:,it),rc,rf,sigref);
    [hfracv0,hfracfv0,kref_fcv0] = get_layer_frac(kref(:,:,:,it),rhoLoc(:,:,:,it),rc,rf,sigref);
    hfrac(:,:,it) = hfrac0; hfracf(:,:,it) = hfracf0; kref_fc(:,:,it) = kref_fc0;
    hfrac_u(:,:,it) = hfracu0; hfracfu(:,:,it) = hfracfu0; kref_fu(:,:,it) = kref_fcu0;
    hfrac_v(:,:,it) = hfracv0; hfracfv(:,:,it) = hfracfv0; kref_fv(:,:,it) = kref_fcv0;
end
%% 2. get the difference of integrals of u and v (dU/dx, dV/dy)

varmom_u = {'uvel_diag'};
varmom_v = {'vvel_diag'};

for ivar = 1:length(varmom_u)
    eval([varmom_u{ivar},'_Int=zeros(nx,ny,nt);']);
end
for ivar = 1:length(varmom_v)
    eval([varmom_v{ivar},'_Int=zeros(nx,ny,nt);']);
end

for it = 1:nt
   disp(['Do integration for u/v, time=  ',num2str(it)])
   for ix = 1:nx
      for iy = 1:ny
          ktmp = kref_fu(ix,iy,it);
          if ktmp~=0
               for ivar = 1:length(varmom_u)
                 eval([varmom_u{ivar},'_Int(ix,iy,it) = sum(',varmom_u{ivar},...
                 '(ix,iy,1:ktmp-1,it).*drF(1,1,1:ktmp-1)) + ',varmom_u{ivar},...
                 '(ix,iy,ktmp,it).*(-hfracfu(ix,iy,it));']);
               end
          end
          ktmp = kref_fv(ix,iy,it);
          if ktmp~=0
               for ivar = 1:length(varmom_v)
                 eval([varmom_v{ivar},'_Int(ix,iy,it) = sum(',varmom_v{ivar},...
                 '(ix,iy,1:ktmp-1,it).*drF(1,1,1:ktmp-1)) + ',varmom_v{ivar},...
                 '(ix,iy,ktmp,it).*(-hfracfv(ix,iy,it));']);
 
               end
          end
      end
   end
end
for ivar = 1:length(varmom_u)
    eval([varmom_u{ivar},'_Int(',varmom_u{ivar},'_Int==0)=nan;']);
end
for ivar = 1:length(varmom_v)
    eval([varmom_v{ivar},'_Int(',varmom_v{ivar},'_Int==0)=nan;']);
end
% get dU/dx and dV/dy
dUdx_Int = zeros(nx,ny);
dVdy_Int = zeros(nx,ny);
% dU/dx
Udy = uvel_diag_Int.*dyc;
dUdx_Int(1:end-1,:) = (Udy(2:end,:)-Udy(1:end-1,:))./raz(2:end,:);
% dV/dy
Vdx = vvel_diag_Int.*dxc;
dVdy_Int(:,1:end-1) = (Vdx(:,2:end)-Vdx(:,1:end-1))./raz(:,2:end);
%%
% % % Um_tot = Um_Diss+Um_Advec+Um_dPH+Um_Ext+Um_PHsur+Um_Impl;
% % % Vm_tot = Vm_Diss+Vm_Advec+Vm_dPH+Vm_Ext+Vm_PHsur+Vm_Impl;
varmom = {'Diss','Adv_a','Cori','PHsur','dPH','Ext','Impl'};
for ivar = 1:length(varmom)
    eval(['Um_',varmom{ivar},'_Int=zeros(nx,ny,nt);']);
    eval(['Vm_',varmom{ivar},'_Int=zeros(nx,ny,nt);']);
end

for it = 1:nt
   disp(['Step 1: Thickness averaging momentum EQ, time=  ',num2str(it)])
   for ix = 1:nx
      for iy = 1:ny
          ktmp = kref_fu(ix,iy,it);
          if ktmp~=0
               for ivar = 1:length(varmom)
                 eval(['Um_',varmom{ivar},'_Int(ix,iy,it) = sum(Um_',varmom{ivar},...
                 '(ix,iy,1:ktmp-1,it).*drF(1,1,1:ktmp-1)) + Um_',varmom{ivar},...
                 '(ix,iy,ktmp,it).*(-hfracfu(ix,iy,it));']);
               end
          end
          ktmp = kref_fv(ix,iy,it);
          if ktmp~=0
               for ivar = 1:length(varmom)
                 eval(['Vm_',varmom{ivar},'_Int(ix,iy,it) = sum(Vm_',varmom{ivar},...
                 '(ix,iy,1:ktmp-1,it).*drF(1,1,1:ktmp-1)) + Vm_',varmom{ivar},...
                 '(ix,iy,ktmp,it).*-(hfracfv(ix,iy,it));']);
               end
          end
      end
   end
end

Um_tot_Int = Um_Diss_Int+Um_Adv_a_Int+Um_Cori_Int+Um_dPH_Int+Um_Ext_Int+Um_PHsur_Int+Um_Impl_Int;
Vm_tot_Int = Vm_Diss_Int+Vm_Adv_a_Int+Vm_Cori_Int+Vm_dPH_Int+Vm_Ext_Int+Vm_PHsur_Int+Vm_Impl_Int;
%% vorticity Calculation
% terms as following
% Vor_Diss, Vor_Adv_a, Vor_Cori, Vor_PH, Vor_Ext, Vor_VVis
disp('Step2: get vorticity of the TWA momentum eq')
Vor_Diss_Int = get_vort_mitgcm_v2(Um_Diss_Int,Vm_Diss_Int,dxc,dyc,raz);
Vor_Adv_a_Int = get_vort_mitgcm_v2(Um_Adv_a_Int,Vm_Adv_a_Int,dxc,dyc,raz);
Vor_Cori_Int = get_vort_mitgcm_v2(Um_Cori_Int,Vm_Cori_Int,dxc,dyc,raz);
% Vor_PH = get_vort_mitgcm_v2(Um_PH_tot,Vm_PH_tot,dxc,dyc,raz);
Vor_PH_00_Int = get_vort_mitgcm_v2(Um_dPH_Int,Vm_dPH_Int,dxc,dyc,raz);
Vor_PHsur_Int = get_vort_mitgcm_v2(Um_PHsur_Int,Vm_PHsur_Int,dxc,dyc,raz);
Vor_Ext_Int = get_vort_mitgcm_v2(Um_Ext_Int,Vm_Ext_Int,dxc,dyc,raz);
Vor_VVis_Int = get_vort_mitgcm_v2(Um_Impl_Int,Vm_Impl_Int,dxc,dyc,raz);

% Vor_tot = Vor_Diss+Vor_Adv_a+Vor_Cori+Vor_PHsur+Vor_PH_00+Vor_Ext+Vor_VVis;
Vor_Diss_Int(Vor_Diss_Int==0)=nan;
Vor_Adv_a_Int(Vor_Adv_a_Int==0)=nan;
Vor_Cori_Int(Vor_Cori_Int==0)=nan;
Vor_PHsur_Int(Vor_PHsur_Int==0)=nan;
Vor_PH_00_Int(Vor_PH_00_Int==0)=nan;
Vor_Ext_Int(Vor_Ext_Int==0)=nan;
Vor_VVis_Int(Vor_VVis_Int==0)=nan;
Vor_tot_Int = Vor_Diss_Int+Vor_Adv_a_Int+Vor_Cori_Int+Vor_PHsur_Int+Vor_PH_00_Int+Vor_Ext_Int+Vor_VVis_Int;
if ichk
    %check u dir.
    amax = 2e-10;
    figure
    subplot(221)
    contourf(Vor_Diss_Int',(-10:1:10)*amax)
    subplot(222)
    contourf(Vor_Adv_a_Int',(-10:1:10)*amax)
    subplot(223)
    contourf(Vor_Cori_Int',(-10:1:10)*amax)
    subplot(224)
    contourf(Vor_PHsur_Int'+Vor_PH_00_Int',(-10:1:10)*amax)
%     subplot(235)
%     contourf(Vor_PH_00_Int',(-10:1:10)*amax)
%     subplot(236)
%     contourf(Vor_Ext_Int',(-10:1:10)*1e-10)
%     subplot(247)
%     contourf(Vor_VVis_Int',(-10:1:10)*1e-10)
%     subplot(248)
%     contourf(Vor_tot_Int',(-10:1:10)*1e-10)
end
    
%% curl of the integrated momentum eqs.
% term1_vort = dZ/dt
% term2_vort = beta*V
% term3_vort = curl(wind)
% term4_vort = curl(bstr)
% term5_vort = f*(dzeta/dt+zeta*dot(b))
% term6_vort = grad(phi_b)xgrad(zeta)
% term7_vort = curl(hvis)
% term8_vort = curl(adv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% term1_vort
% tendency time

% term1_vort
% betav = beta.*vvel_diag_Int;
v_um = Um_Cori_Int./f;
betav = beta.* v_um;
% betav_r = zeros(nx,ny,nt);
% betav_r(:,1:end-1,:) = 0.5*(betav(:,1:end-1,:)+betav(:,2:end,:));
betav_r = 0.5*(betav(2:end,1:end-1)+betav(2:end,2:end));
term2_vort = betav_r;
% term2_vort 
% should be zero
% term3_vort
term3_4_vort = Vor_VVis_Int;
term3_4_vort(isnan(term3_4_vort)) = 0;
% term4_vort
% f(dU/dx+dV/dy)
fdivu = f.*(dUdx_Int+dVdy_Int);
fdivu_psi = 0.25*(fdivu(1:end-1,1:end-1)+fdivu(1:end-1,2:end)+...
    fdivu(2:end,1:end-1)+fdivu(2:end,2:end));
term5_vort = fdivu_psi;
% term5_vort
term6_vort = (Vor_PH_00_Int+Vor_PHsur_Int);
% term6_vort
term7_vort = Vor_Diss_Int;
% term7_vort
term8_vort = Vor_Adv_a_Int;

if ichk
    figure
    plot(-Vor_Cori_Int(:,40))
    hold on
    plot(fdivu_psi(1:end,40)+betav_r(1:end,40))
end
% valid the budget
tot_vort = -term2_vort + term3_4_vort + (-term5_vort) + term6_vort + term7_vort + term8_vort;

figure
plot(tot_vort(:,40),'k','linewidth',2)
hold on

plot(term2_vort(:,40))
plot(term3_4_vort(:,40))
plot(term5_vort(:,40))
plot(term6_vort(:,40))
plot(term7_vort(:,40))
plot(term8_vort(:,40))
%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k0,k0u,k0v] = get_layer_ref1(r0,r0u,r0v,rsig)

nx = size(r0,1);
ny = size(r0,2);

k0 = zeros(nx,ny);
k0u = zeros(nx,ny);
k0v = zeros(nx,ny);
for ix = 1:nx
    for iy = 1:ny
        sigma_t0 = squeeze(r0(ix,iy,:));
        sigma_tu0 = squeeze(r0u(ix,iy,:));
        sigma_tv0 = squeeze(r0v(ix,iy,:));
        ksig0 = find(sigma_t0<=rsig);
        ksigu0 = find(sigma_tu0<=rsig);
        ksigv0 = find(sigma_tv0<=rsig);
        if isempty(ksig0)
            k0(ix,iy) = 0;
        else
            k0(ix,iy) = ksig0(end);
        end
        if isempty(ksigu0)
            k0u(ix,iy) = 0;
        else
            k0u(ix,iy) = ksigu0(end);
        end
        if isempty(ksigv0)
            k0v(ix,iy) = 0;
        else
            k0v(ix,iy) = ksigv0(end);
        end
    end
end

end % function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h1,h2,krefc] = get_layer_frac(kr0,r0,rc,rf,rsig);

nx = size(r0,1);
ny = size(r0,2);
for ix = 1:nx
    for iy = 1:ny
        if kr0(ix,iy)~=0
            htmp1 = rc(kr0(ix,iy));
            htmp2 = rc(kr0(ix,iy)+1);
            rtmp1 = r0(ix,iy,kr0(ix,iy));
            rtmp2 = r0(ix,iy,kr0(ix,iy)+1);
            htmp = htmp1+(rsig-rtmp1)*(htmp2-htmp1)/(rtmp2-rtmp1);
            %              hsec(ix,iy,it) = htmp;
            h1(ix,iy) = -(htmp-htmp1);
            ktmp0 = find(rf>=htmp);
            krefc(ix,iy) = ktmp0(end);
            h2(ix,iy) = htmp-rf(ktmp0(end));
        else
            h1(ix,iy) = 0;
            h2(ix,iy) = 0;
            krefc(ix,iy) = 0;
        end
    end
end

end
