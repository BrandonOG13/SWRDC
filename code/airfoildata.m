function [data] = airfoildata(data)
%This function computes and organizes all airfoil profile, lift and drag data
%   INPUTS: structure DATA
%   OUTPUTS: updated structure DATA

% Aerodynamic Data
file1=importdata(data.airfoildatastring);
foils=reshape(file1,size(file1,1),3,[]);
numfiles=size(foils,3);
aoa=(-180:0.1:180)';

afdata=zeros(length(aoa),3,numfiles);
for i=1:1:numfiles
yi=interp1(foils(:,1,i),foils(:,2:3,i),aoa);
afdata(:,:,i)=[aoa yi];
end
data.afdata=afdata;

%Airfoil Profiles
file2=importdata(data.profiledatastring);
profdata=reshape(file2,size(file2,1),2,[]);

%Profile Data Conditioning
%Normalize the X coordinate data
scale=(max(profdata(:,1))-min(profdata(:,1)));

profdata(:,1)=(profdata(:,1)-min(profdata(:,1)))/...
    scale;
profdata(:,2)=profdata(:,2)/scale;

%Check profdata if CW (LE->TE->LE) and if so, change to CCW (TE->LE->TE)
if (profdata(1,1,1)==0.0 && profdata(end,1,1)==0.0)
ind=find(profdata(:,1,1)==1);
upper=flipdim(profdata(1:ind,:,:),1);
lower=flipdim(profdata(ind:end,:,:),1);
upper(end,:,:)=[];
profdata=cat(1,upper,lower);
end

%Use finer sampled grid
ind2=find(profdata(:,1,1)==0);
xgrid=((cos(0:((2*pi)/300):(2*pi))+1)/2)';
ind3=find(xgrid(:,1,1)==0);

upper=interp1(profdata(1:ind2,1),profdata(1:ind2,2),xgrid(1:ind3),'linear');
lower=interp1(profdata(ind2:end,1),profdata(ind2:end,2),xgrid(ind3:end),'linear');
lower(1,:,:)=[];
profdata=cat(2,xgrid,cat(1,upper,lower));

%Write file for XFOIL analysis
dlmwrite('airfoilXFOIL.txt', profdata, 'precision', '%.6f', ...
    'newline', 'pc')

%Change CCW back to CW (LE->TE->LE)
ind3=find(profdata(:,1,1)==0);
upper=flipdim(profdata(1:ind3,:,:),1);
lower=flipdim(profdata(ind3:end,:,:),1);
upper(end,:,:)=[];
profdata=cat(1,upper,lower);

%Load Data
x=profdata(:,1);
y=profdata(:,2);
limit=find(x(:,1)==1);

%x Coordinates
xu=x(1:limit);
xl=flipud(x(limit:end));

%y Coordinates
yu=y(1:limit,:,:);
yl=flipud(y(limit:end));

%Same number of upper and lower coordinates
numcoords=max(length(xu),length(xl));
xgrid=sort((cos(0:(pi/(numcoords-1)):pi)+1)/2)';
yu=interp1(xu,yu,xgrid);
yl=interp1(xl,yl,xgrid);
xu=xgrid;
xl=xu;

%Export final results
data.limit=limit;
data.yu=yu;
data.yl=yl;
data.xu=xu;
data.xl=xl;

end