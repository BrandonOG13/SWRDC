function [cp,ct,rootflapbend,mass_cons,fq,ts,noise,layouts_dF]=...
    compObjFs(pitch,tsr,twists,chords,data,afdata,nu,rpm_design,rpm_max)
%% Compute Objective Functions
%Run BEM Code
dcp_dct=bem(pitch,tsr,twists,chords,data,afdata,nu);

%Load Results
[cp,ct,rootflapbend,~,~]=loads(dcp_dct,data,nu,rpm_design);

%Calculate AEP instead of CP if specified by user
if data.optAEP==1
    [cp]=aep(cp,data,0);
end

%Calculate Extreme Loads
[pn_worst,pt_worst,tempdata]=iecClass(pitch,twists,chords,afdata,data,nu,rpm_max);

%Structural Calculations at Extreme Loads
add_structure=data.add_structure;
if add_structure==1;
    [mass_cons, dmass, fq, layouts_dF]=struc(chords,twists,pn_worst,pt_worst,tempdata,nu,pitch,data.pvalstruct,rpm_max,add_structure);
elseif data.add_structure==0;
    [mass_cons, dmass, fq, layouts_dF]=struc(chords,twists,pn_worst,pt_worst,tempdata,nu,pitch,data.pvalstruct,rpm_max,add_structure);
end

%Starting Calculation
delt=0.1;
if data.add_starting==1;
    [ts]=starting(data,chords,twists,delt,nu,dmass);
elseif data.add_starting==0;
    ts=NaN;
end

%Noise Calculation
if data.add_noise==1;
    [noise] = noise_calc(chords,twists,nu,data,rpm_design);
elseif data.add_noise==0;
    noise=NaN;
end

end