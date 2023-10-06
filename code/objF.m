function [oc] = objF(data,objcons,zmax,zmin)
%Objectives being considered
zind=data.zind;

%Extract weights
w=data.weights;

%Remove objectives that user has deselected
zmax(zind)=[];
zmin(zind)=[];
objcons(logical([zind 0]))=[];

%Extract objectives
fx=objcons(1:(end-1));

%Scalarization function
singleobj=max(w.*(abs(fx-zmin)./abs(zmax-zmin)));

%Extract constraint violation value
g=objcons(end);

oc=cat(2,singleobj,g);