function [res_vals,chords,twists,thicks] = finalresults(x,data)
%Load required data
nu=data.nu;
afdata=data.afdata;
pitch=data.pitch;
tsr=data.lambda_design;
pval=data.pval;
valBplot=data.valBplot;
printcoords=data.printcoords;
R=data.R;
vo_design=data.vo_design;
rpm_design=(30/pi)*((data.lambda_design*vo_design)/R);%RPM at design lambda and design windspeed
rpm_max=data.rat_rpm;%Max RPM

%Calculate blade geometry
[chords,twists,thicks] = geometry(x',data,nu,pval);

%% Compute Objective Functions
[cp,ct,rootflapbend,mass_cons,fq,ts,noise,layouts_dF]=...
    compObjFs(pitch,tsr,twists,chords,data,afdata,nu,rpm_design,rpm_max);

%% Output Final Blade Results

%Display on command window
cp_ts_m_cons=cat(2,cp,ts,mass_cons,noise);
fprintf('%0.6f % 0.6f % 9.6f % 11.6f % 10.4f % 12.4f % 12.4f % 12.4f\n',cat(2,cp_ts_m_cons,fq))

%Export coordinates
if valBplot==1 || printcoords==1
    [~,~]=exportcoords(R,layouts_dF,valBplot,printcoords,chords,twists,nu);
end

%Calculate results in dimensional quantities
if data.optAEP==0
    power=0.5*data.rho*(vo_design^3)*pi*(R^2)*cp;
    torque=power/(rpm_design*(pi/30));
else
    power=NaN;
    torque=NaN;
end
thrust=0.5*data.rho*(vo_design^2)*pi*(R^2)*ct;

%% Output Results onto Excel
res_vals=cat(2,cp,power,torque,thrust,rootflapbend,ts,mass_cons,fq,noise);

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

%Write inputs in data structure
xlswrite1(saveTextbox,data.numBlades,'Performance Results', 'L12')
xlswrite1(saveTextbox,R,'Performance Results', 'L13')
xlswrite1(saveTextbox,data.hubR,'Performance Results', 'L14')
xlswrite1(saveTextbox,data.rat_rpm,'Performance Results', 'L16')
xlswrite1(saveTextbox,data.J_gen,'Performance Results', 'L17')
xlswrite1(saveTextbox,data.tor_resis,'Performance Results', 'L19')
xlswrite1(saveTextbox,data.rho,'Performance Results', 'L21')
xlswrite1(saveTextbox,data.vis,'Performance Results', 'L22')
xlswrite1(saveTextbox,vo_design,'Performance Results', 'L25')
xlswrite1(saveTextbox,data.lambda_design,'Performance Results', 'L26')
xlswrite1(saveTextbox,rpm_design,'Performance Results', 'L27')

%Write results
xlswrite1(saveTextbox,res_vals,'Performance Results','A2')
xlswrite1(saveTextbox,cat(2,nu*R,chords,twists*(180/pi)),'Optimized Blade Shape','A2')
xlswrite1(saveTextbox,pitch*(180/pi),'Optimized Blade Shape','M8')
if data.cptsr==1
   tsr_rng=linspace(data.tsrmin,data.tsrmax,20);
   [dcp_dct]=bem(pitch,tsr_rng,twists,chords,data,afdata,nu);
   [cp,ct,~,~,~]=loads(dcp_dct,data,nu,rpm_design);
   xlswrite1(saveTextbox, cat(2,tsr_rng',cp,ct),'Optimized Blade Shape','M11')
end

%Write airfoil inputs to Excel file
%collect airfoil data
[h0] = plotairfoildata(data);

%Get name of specified worksheet from workbook
TargetSheet = get(Excel.sheets,'item','Performance Results');

%Paste in the MATLAB figures
print(h0, '-dbitmap'); TargetSheet.Range('A10').PasteSpecial;
delete(h0)%close figure

%Close Excel workbook. 
ExcelWorkbook.Save 
ExcelWorkbook.Close(false)
Excel.Quit; 
delete(Excel); 

end

function [chords,twists] = exportcoords(R,layouts_d,valBplot,printcoords,chords,twists,nu)
%Exports blade geometry coordinates for input into CAD
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com

%% Plot/Save Results

layouts_d_ori=layouts_d;
layouts_d(1,:,:)=[];%Used for cap since LE is susceptible to intersection
numlayers=2;
thDim=size(layouts_d,3);
a=(0:numlayers:(thDim-numlayers));

red=1; %Increase value to reduce number of points in 3D blade plot
limit=size(layouts_d,1);
limit_ori=size(layouts_d_ori,1);
coarse=(1:red:limit)';
coarse_ori=(1:red:limit_ori)';
lencoarse=length(coarse);
lencoarse_ori=length(coarse_ori);
sq_pu=zeros(lencoarse_ori,numlayers);
sq_pl=zeros(lencoarse_ori,numlayers);
sq_u=zeros(lencoarse_ori,numlayers);
sq_l=zeros(lencoarse_ori,numlayers);
spar_xu=zeros(lencoarse,numlayers);
spar_xl=zeros(lencoarse,numlayers);
spar_yu=zeros(lencoarse,numlayers);
spar_yl=zeros(lencoarse,numlayers);

for i=1:length(a)
    j=1;
    xp=layouts_d(:,1,j+a(i));
    comp=numlayers;
    xspar=layouts_d(:,1,a(i)+comp);
    xsparco=interp1((1:limit)',xspar,coarse);
    yuspar=interp1((1:limit)',layouts_d(:,2,a(i)+comp),coarse);
    ylspar=interp1((1:limit)',layouts_d(:,3,a(i)+comp),coarse);
    xspar_u=(xsparco*cos(twists(i)))-(yuspar*sin(twists(i)));
    yspar_u=(xsparco*sin(twists(i)))+(yuspar*cos(twists(i)));
    xspar_l=(xsparco*cos(twists(i)))-(ylspar*sin(twists(i)));
    yspar_l=(xsparco*sin(twists(i)))+(ylspar*cos(twists(i)));
    spar_xu(:,i)=xspar_u;
    spar_xl(:,i)=xspar_l;
    spar_yu(:,i)=yspar_u;
    spar_yl(:,i)=yspar_l;
    
    xp_ori=layouts_d_ori(:,1,j+a(i));
    xpc=interp1((1:limit_ori)',xp_ori,coarse_ori);
    yuc=interp1((1:limit_ori)',layouts_d_ori(:,2,j+a(i)),coarse_ori);
    ylc=interp1((1:limit_ori)',layouts_d_ori(:,3,j+a(i)),coarse_ori);
    
    xru=(xpc*cos(twists(i)))-(yuc*sin(twists(i)));
    yru=(xpc*sin(twists(i)))+(yuc*cos(twists(i)));
    xrl=(xpc*cos(twists(i)))-(ylc*sin(twists(i)));
    yrl=(xpc*sin(twists(i)))+(ylc*cos(twists(i)));
    
    sq_pu(:,i)=xru;
    sq_pl(:,i)=xrl;
    sq_u(:,i)=yru;
    sq_l(:,i)=yrl;
end
lenR=repmat(nu'*R,lencoarse,1);
lenR_ori=repmat(nu'*R,lencoarse_ori,1);

%% Display 3D Blade Shape
if valBplot==1
    %If figure was already open, close it
    b3d=figure(10);
    close(b3d)
    
    %Generate the figure
    b3d=figure(10);
    hold on
    surf(lenR,spar_xu,spar_yu,'EdgeColor','none','FaceColor',[0 0 1])
    surf(lenR,spar_xl,spar_yl,'EdgeColor','none','FaceColor',[0 0 1])
    
    surf(lenR_ori,sq_pu,sq_u,'EdgeColor','none','FaceAlpha',0.8)
    surf(lenR_ori,sq_pl,sq_l,'EdgeColor','none','FaceAlpha',0.8)
    colormap gray
    
    set(gca,'position',[0 0 1 1]);
    set(gca,'visible','off')
    hold off
    set(b3d, 'units', 'centimeters', 'pos', [1 1 12 8])
    axis image;
    view([-55,15])
end

%% Export for Rapid Prototyping
if printcoords==1
    %x and y cap
    xC=cat(1,spar_xu,flipdim(spar_xl,1),spar_xu(1,:));
    yC=cat(1,spar_yu,flipdim(spar_yl,1),spar_yu(1,:));
    
    xC(limit,:)=[];
    yC(limit,:)=[];
    
    %x and y profile
    xP=cat(1,sq_pu,flipdim(sq_pl,1));
    yP=cat(1,sq_u,flipdim(sq_l,1));
    
    xP(limit_ori,:)=[];
    yP(limit_ori,:)=[];
    
    %z
    limit2=size(xC,1);
    limitP=size(xP,1);
    rsB=(nu-min(nu))'*R;
    
    z=reshape(repmat(rsB,limit2,1),limit2,1,[]);
    zP=reshape(repmat(rsB,limitP,1),limitP,1,[]);
    
    %Combine
    xyzC=cat(2,reshape(xC,limit2,1,[]),reshape(yC,limit2,1,[]),z);
    xyzP=cat(2,reshape(xP,limitP,1,[]),reshape(yP,limitP,1,[]),zP);
    
    %Write to text files
    for i=1:size(xyzC,3)
        str2=cat(2,'exportedcoords\cap_coords',num2str(i),'.txt');
        dlmwrite(str2,xyzC(:,:,i), 'delimiter', '\t','precision', 6)
    end
    for i=1:size(xyzP,3)
        str1=cat(2,'exportedcoords\pro_coords',num2str(i),'.txt');
        dlmwrite(str1,xyzP(:,:,i), 'delimiter', '\t','precision', 6)
    end
    dlmwrite('exportedcoords\final_blade_shape.txt',cat(2,(nu*R),chords,twists*(180/pi)),'\t')
end

end