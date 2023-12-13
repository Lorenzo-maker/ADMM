function [Fx0_fun, Fy0_fun, Fx_tyre_fun, Fy_tyre_fun, Fxlim_fun, Fylim_fun, Gxa_fun, Gyk_fun, mf_fun] = MF_Tyre(tir, options)
    arguments
        tir
        options.shift = 1;
        options.curvature = 1;
        options.sym = 'SX';
    end
    curv = options.curvature;
    shift = options.shift;
    import casadi.*
    tol = 1e-4;
    Fz = SX.sym('Fz'); % tyre vertical load 
    dfz = (Fz - tir.FNOMIN)/tir.FNOMIN;
    gamma = SX.sym('gamma');
    
    % Pag. 176/3-th edition
    % Longitudinal Force 
    % Vsx = Vx - omega*r (Vx = vel. centro ruota)
    % Vr = Vx - Vsx = omega*r    
    k = SX.sym('k'); % longitudinal slip (= - Vsx/|Vx|)
    Shx = shift*(tir.PHX1 + tir.PHX2*dfz);
    kx = k + Shx;
    Cx = tir.PCX1;
    mu_x = (tir.PDX1 + tir.PDX2*dfz)*(1 - tir.PDX3*gamma^2); % No data for V0 and Vs (i.e. Vs/V0 = 1)
    Dx = mu_x*Fz;
    Ex = (tir.PEX1 + tir.PEX2*dfz + tir.PEX3*dfz^2)*(1 - curv*tir.PEX4*if_else_smooth(kx, 0, 1, -1, 'C', 1e4));
    Kxk = Fz*(tir.PKX1 + tir.PKX2*dfz)*exp(tir.PKX3*dfz);
    Bx = Kxk/(Cx*Dx + tol);
    Svx = shift*(Fz*(tir.PVX1 + tir.PVX2*dfz)); % approx (|Vx|/(eps_V + |Vx|) = 1)
    Fx0 = -Dx*sin(Cx*atan(Bx*kx - Ex*(Bx*kx - atan(Bx*kx)))) + Svx;
    Fx0_fun = Function('Fx0', {k, Fz, gamma}, {Fx0});
    Fxmin = -Dx + shift*Svx;
    Fxmax = Dx + shift*Svx;    
    Fxlim_fun = Function('Fxlim', {Fz, gamma}, {Fxmin, Fxmax});
    
    % Lateral Force 
    % Vsy = Vy (vel. laterale centro ruota)
    % tg(alfa) = - Vy/Vx
    alpha = SX.sym('alpha');
    Kya = tir.PKY1*tir.FNOMIN*sin(tir.PKY4*atan(Fz/(tir.FNOMIN*(tir.PKY2 + tir.PKY5*gamma^2))))*(1 - tir.PKY3*if_else_smooth(gamma,0,gamma,-gamma, 'C', 1000));
    Svyg = shift*(Fz*(tir.PVY3 + tir.PVY4*dfz)*gamma);
    Svy = shift*(Fz*(tir.PVY1 + tir.PVY2*dfz) + Svyg);
    Ky0g = Fz*(tir.PKY6 + tir.PKY7*dfz);
    Shy = (tir.PHY1 + tir.PHY2*dfz) + (Ky0g*gamma - Svyg)/(Kya + tol);
    alphay = alpha + shift*Shy;
    Cy = tir.PCY1;
    mu_y = (tir.PDY1 + tir.PDY2*dfz)*(1 - tir.PDY3*gamma^2);
    Dy = mu_y*Fz;
    Ey = (tir.PEY1 + tir.PEY2*dfz)*(1 + tir.PEY5*gamma^2 - curv*(tir.PEY3 + tir.PEY4*gamma)*if_else_smooth(alphay, 0, 1, -1, 'C', 1e4));
    By = Kya/(Cy*Dy + tol);
    Fy0 = -Dy*sin(Cy*atan(By*alphay - Ey*(By*alphay - atan(By*alphay)))) + Svy;
    Fy0_fun = Function('Fy0', {alpha, Fz, gamma}, {Fy0});
    Fymin = -Dy + Svy;
    Fymax = Dy + Svy;
    Fylim_fun = Function('Fymax', {Fz, gamma}, {Fymin, Fymax});
    
    % Longitudinal Force 
    Bxa = (tir.RBX1 + tir.RBX3*gamma^2)*cos(atan(tir.RBX2*k));
    Cxa = tir.RCX1;
    Exa = tir.REX1 + tir.REX2*dfz;
    Shxa = shift*(tir.RHX1);
    alphas = (alpha + Shxa);
    Gxa0 = cos(Cxa*atan(Bxa*Shxa - Exa*(Bxa*Shxa - atan(Bxa*Shxa))));
    Gxa = cos(Cxa*atan(Bxa*alphas - Exa*(Bxa*alphas - atan(Bxa*alphas))))/Gxa0;
    Gxa_fun = Function('Gxa_tyre', {k, alpha, Fz, gamma}, {Gxa});
    Fx_tyre = Gxa*Fx0_fun(k, Fz, gamma);
    Fx_tyre_fun = Function('Fx_tyre', {k, alpha, Fz, gamma}, {Fx_tyre});
    
    % Longitudinal Force 
    Shyk = shift*(tir.RHY1 + tir.RHY2*dfz);
    ks = k + Shyk;
    Byk = (tir.RBY1 + tir.RBY4*gamma^2)*cos(atan(tir.RBY2*(alpha - tir.RBY3)));
    Cyk = tir.RCY1;
    Eyk = tir.REY1 + tir.REY2*dfz;
    Dvyk = mu_y*Fz*(tir.RVY1 + tir.RVY2*dfz + tir.RVY3*gamma)*cos(atan(tir.RVY4*alpha));
    Svyk = shift*(Dvyk*sin(tir.RVY5*atan(tir.RVY6*k)));
    Gyk0 = cos(Cyk*atan(Byk*Shyk - Eyk*(Byk*Shyk - atan(Byk*Shyk))));
    Gyk = cos(Cyk*atan(Byk*ks - Eyk*(Byk*ks - atan(Byk*ks))))/Gyk0;
    Gyk_fun = Function('Gyk_fun', {alpha, k, Fz, gamma}, {Gyk});
    Fy_tyre = Gyk*Fy0_fun(alpha, Fz, gamma) + Svyk;
    Fy_tyre_fun = Function('Fy_tyre', {alpha, k, Fz, gamma}, {Fy_tyre});
    mf_fun = casadi.Function('mf',{k,alpha,Fz,gamma},{Fx_tyre,Fy_tyre},struct('cse',true));


    
end
