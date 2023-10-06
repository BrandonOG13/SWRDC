function [data] = collect_data(handles,testflag)
%This function collects all the data specified on the GUI by the user

%Design data
data.numBlades=str2num(get(handles.numBlades,'String'));
data.pitch=str2num(get(handles.pitch,'String'));
data.vo_design=str2num(get(handles.vo_design,'String'));
data.R=str2num(get(handles.R,'String'));
data.hubR=str2num(get(handles.hubR,'String'));
data.rat_rpm=str2num(get(handles.rpm,'String'));
data.lambda_design=str2num(get(handles.lambda_design,'String'));
if ((data.lambda_design*data.vo_design)/data.R)>(data.rat_rpm*(pi/30))
    disp('ABORT: Design RPM exceeds the maximum allowable RPM of generator')
    data=NaN;
    return
end
data.maxpower=str2num(get(handles.maxpower,'String'));

%Starting
data.u_start=str2num(get(handles.u_start,'String'));
data.lambda_start=str2num(get(handles.lambda_start,'String'));
data.J_gen=str2num(get(handles.Jgen,'String'));
data.tor_resis=str2num(get(handles.tor_resis,'String'));
data.add_starting=get(handles.add_starting,'Value');

%Wind data
data.rho=str2num(get(handles.rho,'String'));
data.vis=str2num(get(handles.vis,'String'));

%AEP data
data.optAEP=get(handles.optAEP,'Value');
data.winddata=get(handles.winddata,'String');
data.mws=str2num(get(handles.mws,'String'));
data.probtype=get(handles.probtype,'Value');
data.vcutin=str2num(get(handles.vcutin,'String'));
data.vcutout=str2num(get(handles.vcutout,'String'));
data.k_weibull=str2num(get(handles.k_weibull,'String'));
data.A_weibull=str2num(get(handles.A_weibull,'String'));

%BEM
data.iter=str2num(get(handles.iter,'String'));
data.belements=str2num(get(handles.belements,'String'));
data.relax=str2num(get(handles.relax,'String'));

%Plotting
data.pval=get(handles.pval,'Value');
data.pvalstruct=get(handles.pvalstruct,'Value');
data.valBplot=get(handles.valBplot,'Value');
data.optprogress=get(handles.optprogress,'Value');
data.printcoords=get(handles.printcoords,'Value');
data.tsrmax=str2num(get(handles.tsrmax,'String'));
data.tsrmin=str2num(get(handles.tsrmin,'String'));
data.cptsr=get(handles.cptsr,'Value');

%Structural data
data.L=str2num(get(handles.L,'String'));
data.U=str2num(get(handles.U,'String'));
data.E=str2num(get(handles.E,'String'));
data.rh=str2num(get(handles.rh,'String'));
data.ep_failure=str2num(get(handles.ep_failure,'String'));
data.materialSF=str2num(get(handles.materialSF,'String'));
data.MinsparTh=str2num(get(handles.MinsparTh,'String'));
data.add_structure=get(handles.add_structure,'Value')-1;
data.perform_structure=get(handles.perform_structure,'Value');

%Blade parameterization data
data.N=str2num(get(handles.N,'String'));
data.uptwist=str2num(get(handles.uptwist,'String'));
data.dotwist=str2num(get(handles.dotwist,'String'));
data.upchord=str2num(get(handles.upchord,'String'));
data.dochord=str2num(get(handles.dochord,'String'));
data.spacing=get(handles.spacing,'Value');
data.spacingBE=get(handles.spacingBE,'Value');

%Airfoil data
data.airfoildatastring=(get(handles.airfoildata,'String'));
data.profiledatastring=(get(handles.airfoilprofile,'String'));
data.re_array=str2num(get(handles.re_array,'String'));
data.thick=str2num(get(handles.thick,'String'));

%Optimization Data
w1=str2num(get(handles.w1,'String'));
w2=str2num(get(handles.w2,'String'));
w3=str2num(get(handles.w3,'String'));
w4=str2num(get(handles.w4,'String'));
data.weights=cat(2,w1,w2,w3,w4);

%Genetic Algorithm Data
data.pop=str2num(get(handles.pop,'String'));
if rem(data.pop/2,2)~=0
    disp('ABORT: Population size divided by 2 must be an even number')
    data=NaN;
    return
end
data.mingen=str2num(get(handles.mingen,'String'));
data.nm=str2num(get(handles.nm,'String'));
data.pm=str2num(get(handles.pm,'String'));
data.val=get(handles.val,'Value');
data.nc=str2num(get(handles.nc,'String'));
data.pc=str2num(get(handles.pc,'String'));
data.flag3=get(handles.flag3,'Value');

%Noise data
data.zo=get(handles.zo,'Value');
data.Ratio=str2num(get(handles.Ratio,'String'));
data.a_top=str2num(get(handles.a_top,'String'));
data.a_ground=str2num(get(handles.a_ground,'String'));
data.Htower=str2num(get(handles.Htower,'String'));
data.shaftlength=str2num(get(handles.shaftlength,'String'));
data.Gamma=str2num(get(handles.Gamma,'String'));
data.estimate=get(handles.estimate,'Value');
data.Tscale=str2num(get(handles.Tscale,'String'));
data.Tinten=str2num(get(handles.Tinten,'String'));
data.r0=str2num(get(handles.r0,'String'));
data.round=get(handles.round,'Value');
if data.round==2; data.round=0; end
data.hblunt=str2num(get(handles.hblunt,'String'));
data.TEangle=str2num(get(handles.TEangle,'String'));
data.h0=str2num(get(handles.h0,'String'));
data.PSI=str2num(get(handles.PSI,'String'));
data.add_noise=get(handles.add_noise,'Value');
data.noisedatabase=(get(handles.noisedatabase,'String'));
data.Rref=str2num(get(handles.Rref,'String'));
data.Aref=str2num(get(handles.Aref,'String'));

%XFOIL data
data.ITER=str2num(get(handles.ITER,'String'));
data.Ncrit=get(handles.Ncrit,'Value');
data.xfoilexecute=get(handles.xfoilexecute,'Value');
if data.add_noise==1
    if data.xfoilexecute==0
        try    temp=importdata(data.noisedatabase);    catch
            disp('ABORT: Noise database file not found.');    data=NaN;
            return
        end
    end
end

%Load data
data.loadcase=get(handles.loadcases, 'Value');
data.rotoroverspeed=str2num(get(handles.rotoroverspeed,'String'));
data.Ve50=str2num(get(handles.Ve50,'String'));
data.loadSF=str2num(get(handles.loadSF,'String'));

%Baseline Blade data
data.turbdata=(get(handles.turbdata,'String'));
data.runBaseline=get(handles.runBaseline,'Value');
if data.runBaseline==1
    try    temp=importdata(data.turbdata);    catch
        disp('ABORT: Baseline blade data file not found.');    data=NaN;
        return
    end
end

%Parallel Computing Toolbox - MATLAB - MathWorks
data.parallel=get(handles.parallel,'Value');

%Save image of GUI 
set(gcf,'PaperPositionMode','auto')
saveas(gcf,'GUI_settings\GUI_Input.bmp','bmp');

%Filename and folder specifying results file
data.saveTextbox=get(handles.saveTextbox,'String');

%Run the setup file
data=setup(data,1);
end

