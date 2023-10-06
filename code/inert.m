function [stiffs_mass_cg]=inert(layouts,data,numlayers)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
xc=squeeze(layouts(:,1,:));
yu=squeeze(layouts(:,2,:));
yl=squeeze(layouts(:,3,:));

x1=flipdim((cat(1,xc,flipdim(xc,1))),1);
x0=cat(1,x1(end,:),x1(1:(end-1),:));

y1=flipdim((cat(1,yu,flipdim(yl,1))),1);
y0=cat(1,y1(end,:),y1(1:(end-1),:));

area=0.5*((x0.*y1)-(x1.*y0));

qx=(1/6)*(x0-x1).*((y0.^2)+(y0.*y1)+(y1.^2));
qy=(1/6)*(y1-y0).*((x0.^2)+(x0.*x1)+(x1.^2));

ixx=(1/12)*(x0-x1).*(y0+y1).*((y0.^2)+(y1.^2));
iyy=(1/12)*(y1-y0).*(x0+x1).*((x0.^2)+(x1.^2));
ixy=(1/24)*((x0.*y1)-(x1.*y0)).*((2*x0.*y0)+(x0.*y1)+(x1.*y0)+(2*x1.*y1));

moms_out = squeeze(sum(cat(3,area,qx,qy,ixx,iyy,ixy),1));

belements=size(layouts,3)/numlayers;
props=reshape(moms_out,[numlayers belements 6]);

%shell
eq=(props(1,:,:)-props(2,:,:));

trunc_props_E=cat(1,eq*data.E);
trunc_props_rh=cat(1,eq*data.rh);

%Stiffness about Reference Axes with origin set at Leading Edge
EI_flap=sum(trunc_props_E(:,:,4),1);
EI_edge=sum(trunc_props_E(:,:,5),1);
EI_prod=sum(trunc_props_E(:,:,6),1);
LE_theta=(0.5*atan((2*EI_prod)./(EI_edge-EI_flap)));
LE_theta=round(LE_theta*(10^5))/(10^5);
LE_theta(isnan(LE_theta))=0;

%Center of Mass
mass=sum(trunc_props_rh(:,:,1),1);
mQXsec=sum(trunc_props_rh(:,:,2),1);
mQYsec=sum(trunc_props_rh(:,:,3),1);
mXEvec=(mQYsec./mass);
mYEvec=(mQXsec./mass);

%Point of Elasticity
EAsec=sum(trunc_props_E(:,:,1),1);
EQXsec=sum(trunc_props_E(:,:,2),1);
EQYsec=sum(trunc_props_E(:,:,3),1);

eXEvec=(EQYsec./EAsec);
eYEvec=(EQXsec./EAsec);

EIxSUM=EI_flap-((eYEvec.^2).*EAsec);
EIySUM=EI_edge-((eXEvec.^2).*EAsec);
EIxySUM=EI_prod-((eXEvec.*eYEvec).*EAsec);

theta=(0.5*atan((2*EIxySUM)./(EIySUM-EIxSUM)));
theta(isnan(theta))=0;

EI1=((0.5*(EIxSUM+EIySUM))+(0.5*(EIxSUM-EIySUM).*cos(2*theta))-(EIxySUM.*sin(2*theta)));
EI2=((0.5*(EIxSUM+EIySUM))-(0.5*(EIxSUM-EIySUM).*cos(2*theta))+(EIxySUM.*sin(2*theta)));

stiffs_mass_cg=cat(2,EI1',EI2',theta',mass',mXEvec',mYEvec',eXEvec',eYEvec',EI_flap',EI_edge',EAsec',LE_theta',EIxSUM',EIySUM',EIxySUM');