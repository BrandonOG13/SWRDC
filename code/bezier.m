function [bedim] = bezier(ctrlpnts,N,nu,cp_radii)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
%Compute Bernstein Polynomials
t1=(nu-(min(nu)));
t=(t1/max(t1));
length_t=length(t);
i=(0:1:N);
imat=repmat(i,length_t,1);
bino=factorial(N)./(factorial(imat).*factorial(N-imat));
B=zeros(length_t,(N+1));
for j=1:1:(N+1)
    B(:,j)=(t.^i(j)).*((1-t).^(N-i(j)));
end
B=bino.*B;

length_t=size(B,1);

values=sum((B.*repmat(ctrlpnts',length_t,1)),2);
radl=sum((B.*repmat(cp_radii,length_t,1)),2);      %Radii Locations

bedim=interp1(radl,values,nu);
end