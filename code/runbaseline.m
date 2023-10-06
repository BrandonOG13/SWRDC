function [cp_ts_m_cons] = runbaseline(data)
%% Load Turbine Data
turbdata=importdata(data.turbdata);
turbdata=turbdata.data;
afdata=data.afdata;

nu_base=turbdata(:,1)/data.R;
chords=turbdata(:,2);
tip=(turbdata(end,3));
%Set tip twist to zero and adjust distribution accordingly
if tip<=0
    twists=(turbdata(:,3)+abs(tip))*(pi/180);
elseif tip>0
    twists=(turbdata(:,3)-abs(tip))*(pi/180);
end
pitch=turbdata(end,3)*(pi/180);

R=data.R;
vo_design=data.vo_design;
tsr=data.lambda_design;
lambda_design=data.lambda_design;

rpm_design=(30/pi)*((lambda_design*vo_design)/R);%RPM at design lambda and design windspeed
rpm_max=data.rat_rpm;%Max RPM

% Display geometry plot if user selected
if data.pval==1
    figure(1)
    subplot(2,2,1);
    plot(nu_base,(twists*(180/pi)),'-b.')
    xlabel('r/R(dimensionless)');   ylabel('Twist (degrees)');
    legendoptions1=legend('Twist');
    set(legendoptions1,'Location','NorthEast');
    title('Twist vs. r/R')
    
    subplot(2,2,2);
    plot(nu_base,chords,'-r.')
    xlabel('r/R(dimensionless)');   ylabel('Chord (m)');
    legendoptions1=legend('Chord');
    set(legendoptions1,'Location','NorthEast');
    title('Chord vs. r/R')
    
    subplot(2,2,3);
    plot(nu_base,(data.thick.*ones(size(nu_base))),'-g.')
    xlabel('r/R(dimensionless)');   ylabel('Thickness/Chord (%)');
    legendoptions1=legend('Thickness');
    set(legendoptions1,'Location','NorthEast');
    title('% Thickness vs. r/R')
    
    subplot(2,2,4);
    plot(nu_base,chords.*data.thick,'-k.')
    xlabel('r/R(dimensionless)');   ylabel('Thickness (m)');
    legendoptions1=legend('Thickness');
    set(legendoptions1,'Location','NorthEast');
    title('Dimensional Thickness vs. r/R')
end

%% Display Header for Baseline Run
disp(' ')
fprintf('Loading and calculating baseline design...\n')
disp(' ')
disp('Cp / AEP  Starting  Mass (kg)  C.1: Strain  C.2: Freq.L  C.3: Freq.U  Noise (dB)  1stFlapEigenFreq (rad/s)')

%% Compute Objective Functions
[cp,ct,rootflapbend,mass_cons,fq,ts,noise,~]=...
    compObjFs(pitch,tsr,twists,chords,data,afdata,...
    nu_base,rpm_design,rpm_max);

%% Output Results
%Display result on command window
cp_ts_m_cons=cat(2,cp,ts,mass_cons,noise);
fprintf('%0.6f % 0.6f % 9.6f % 11.6f % 10.4f % 12.4f % 12.4f % 12.4f\n',cat(2,cp_ts_m_cons,fq))
disp(' ')

power=0.5*data.rho*(vo_design^3)*pi*(R^2)*cp;
torque=power/(rpm_design*(pi/30));
thrust=0.5*data.rho*(vo_design^2)*pi*(R^2)*ct;

%Output Results
res_vals=cat(2,cp,power,torque,thrust,rootflapbend,ts,mass_cons,fq,noise);

%% Write on Excel File
saveTextbox=data.saveTextbox;
%Check if Excel file is open before writing
[fid, ~] = fopen(saveTextbox,'a');
if fid==-1
    disp('ABORT: PLEASE MAKE SURE EXCEL IS NOT RUNNING IN THE BACKGROUND.')
    %disp('Excel results file is already open. Closing file.')
    %wbkname = 'results.xls';
    %h = actxGetRunningServer('Excel.Application');
    %h.WorkBooks.Item(wbkname).Close;
    return
else
    fclose(fid);
end

%Open activex server and check if file already exists,
%create it if it doesn't
Excel = actxserver ('Excel.Application');
File=saveTextbox;
if ~exist(File,'file')
    ExcelWorkbook = Excel.workbooks.Add;
    ExcelWorkbook.SaveAs(File,1);
    ExcelWorkbook.Close(false);
end
ExcelWorkbook = Excel.workbooks.Open(File);

xlswrite1(saveTextbox, res_vals, 'Performance Results', 'A3')
xlswrite1(saveTextbox, cat(2,nu_base*R,chords,twists*(180/pi)),'Baseline Blade Shape','A2')
xlswrite1(saveTextbox, pitch*(180/pi),'Baseline Blade Shape','M8')
if data.cptsr==1
    tsr_rng=linspace(data.tsrmin,data.tsrmax,20);
    [dcp_dct]=bem(pitch,tsr_rng,twists,chords,data,afdata,nu_base);
    [cp,ct,~,~,~]=loads(dcp_dct,data,nu_base,rpm_design);
    xlswrite1(saveTextbox, cat(2,tsr_rng',cp,ct),'Baseline Blade Shape','M11')
end

%Close Excel workbook.
ExcelWorkbook.Save
ExcelWorkbook.Close(false)
Excel.Quit;
delete(Excel);

end
