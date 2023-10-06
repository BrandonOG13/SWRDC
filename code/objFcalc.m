%% AERODYNAMIC, STARTING, STRUCTURAL OPTIMIZATION
function [fx_cons] = objFcalc(chromo,data)
% credit for this source code goes to:
% Matias Sessarego, M.Sc., msessare@ucalgary.ca, msessare@gmail.com
[cp_ts_m_cons, ~]=aero_start_struct(chromo',data);

%Objective functions
f1=cp_ts_m_cons(1);
f2=cp_ts_m_cons(2);
f3=cp_ts_m_cons(3);
f4=cp_ts_m_cons(7);

%Normalize Constraint Violation Values
g=zeros(4,1);
g(1)=cp_ts_m_cons(4);
g(2)=cp_ts_m_cons(5);
g(3)=cp_ts_m_cons(6);

%Check penalty on starting
t_max=data.t_max;
if (f2==((5.5*t_max)+1))||...
        (f2==((4.5*t_max)+1))||...
        (f2==((3.5*t_max)+1))
    g(4)=f2;
end
%Check penalty on power
if isinf(f1)
    g(4)=f1;
end

%Set Feasible Constraint Values to Zero
g(g<=0)=0;
%--------------------------------------

g=sum(abs(g));

fx_cons=cat(2,f1,f2,f3,f4,g);