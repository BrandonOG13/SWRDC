function time_to_start = starting(data,chord,twist,delt,nu,dmass)

u_start=data.u_start;
t_max=data.t_max;
lambda_start=data.lambda_start;
tor_resis=data.tor_resis;

% Determine Inertia
J = data.numBlades*trapz(nu*data.R,(dmass.*((nu*data.R).^2))) + data.J_gen;

% Determine the two constants (c1, c2) in the equation for
% d(lambda)/dt

tw_rad = twist; % Twist must be in radians for starting calcs

% constant terms in front of equation 4.8
c1 = data.rho*data.numBlades*u_start^2*data.R^3;
% constant terms in front of equation 4.9
c2 = data.R/(J*u_start);
lambda = 0.0; time = 0.0;
Fn_3=0.0; Fn_2=0.0; Fn_1=0.0;

while (1) % Use Adams Moulton for integration
    fn_0 = deriv(chord, tw_rad, lambda, data, nu);
    Fn_0 = c2*(c1*fn_0 - tor_resis);
    lam_pred=lambda + delt*(55*Fn_0-59*Fn_1+37*Fn_2-9*Fn_3)/24;
    fp_1 = deriv(chord, tw_rad, lam_pred, data, nu);
    Fp_1 = c2*(c1*fp_1 - data.tor_resis);
    del_l = delt*(9*Fp_1+19*Fn_0-5*Fn_1+Fn_2)/24;
    if (lambda + del_l > lambda_start)
        time_to_start = time + delt*(lambda_start - lambda)/del_l;
        return
    end
    lambda = lambda + del_l;
    time = time + delt  ;
    
    if (time > t_max) % Penalise slow blades
        time_to_start = 5.5*t_max + 1.0 ;
        return
    end
    if (fp_1 < -1)
        time_to_start = 4.5*t_max + 1.0;
        return
    end
    if (fp_1 < 0.0 | lambda < 0)
        time_to_start = 3.5*t_max + 1.0;
        return
    end
    
    Fn_3=Fn_2;
    Fn_2=Fn_1;
    Fn_1=Fn_0;
end
end %start_calc

function dlamdt = deriv(chord, twist, lambda, data, nu)
rad=nu;
lamr = lambda.*rad;
lamr2=lamr.*lamr;
tmp = sqrt(1 + lamr2).*sin(twist).*(cos(twist) - lamr.*sin(twist)).*...
		chord.*rad;
dlamdt=trapz(nu*data.R,tmp);
end