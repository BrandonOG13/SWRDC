function [cp,ct,rootflapbend,pn,pt]=loads(dcp_dct,data,nu,rpm)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
%% Power, Pitch, RPM, Thrust Curves
%Determine optimum pitch & tsr for variable speed fixed pitch region
R=data.R;
rho=data.rho;
numBlades=data.numBlades;
pitch=data.pitch;
omega=rpm*(pi/30);
dcp=dcp_dct(:,:,1);
dct=dcp_dct(:,:,2);
cp=trapz(nu,dcp,2);
ct=trapz(nu,dct,2);

if any(cp>(16/27));    cp(cp>(16/27))=0;     end
if any(cp<0);          cp(cp<0)=0;     end

pn=((dct'*rho*(data.vo_design^2)*pi*R)/(2*numBlades));
pt=(dcp'*rho*(data.vo_design^3)*pi)./(2*omega*repmat(nu',size(dcp,1),1)*numBlades)';

nurep=repmat(nu,1,size(pn,2));
%M - Root Flap Bending Moment
rootflapbend=trapz((nu*R),((nurep*R).*pn));

end