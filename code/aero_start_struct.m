function [cp_ts_m_cons, layouts_dF]=aero_start_struct(x,data)
% credit for this source code goes to:
% Matias Sessarego

%% Load Turbine Data
nu=data.nu;
afdata=data.afdata;
pitch=data.pitch;
tsr=data.lambda_design;
pval=data.pval;

rpm_design=(30/pi)*((data.lambda_design*data.vo_design)/data.R);%RPM at design lambda and design windspeed
rpm_max=data.rat_rpm;%Max RPM

%Calculate blade geometry
[chords,twists,~] = geometry(x,data,nu,pval);

%% Compute Objective Functions
[cp,~,~,mass_cons,fq,ts,noise,layouts_dF]=...
    compObjFs(pitch,tsr,twists,chords,data,afdata,nu,rpm_design,rpm_max);

%% Output Results
cp_ts_m_cons=cat(2,1/cp,ts,mass_cons,noise);

%Display result on command window
cp_ts_m_cons_disp=cp_ts_m_cons;
cp_ts_m_cons_disp(1)=(1/cp_ts_m_cons_disp(1));
fprintf('%0.6f % 0.6f % 9.6f % 11.6f % 10.4f % 12.4f % 12.4f % 12.4f\n',cat(2,cp_ts_m_cons_disp,fq))