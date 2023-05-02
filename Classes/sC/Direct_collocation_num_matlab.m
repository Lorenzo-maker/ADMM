function x_next = Direct_collocation_num(x, u, xdot, dt, options)


arguments
   x
   u
   xdot
   dt
   options.degree = 2;
end
%% Matrix for direct collocation
    import casadi.*

    % Degree of interpolating polynomial
    d = options.degree;

    % Get collocation points
    tau_root = [0 collocation_points(d, 'radau')];

    % Coefficients of the collocation equation
    C = zeros(d+1,d+1);

    % Coefficients of the continuity equation
    D = zeros(d+1, 1);

    % Coefficients of the quadrature function
    B = zeros(d+1, 1);

    % Construct polynomial basis
    for j=1:d+1
      % Construct Lagrange polynomials to get the polynomial basis at the collocation point
      coeff = 1;
      for r=1:d+1
        if r ~= j
          coeff = conv(coeff, [1, -tau_root(r)]);
          coeff = coeff / (tau_root(j)-tau_root(r));
        end
      end
      % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D(j) = polyval(coeff, 1.0);

      % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      pder = polyder(coeff);
      for r=1:d+1
        C(j,r) = polyval(pder, tau_root(r));
      end

      % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      pint = polyint(coeff);
      B(j) = polyval(pint, 1.0);
    end

h = dt;
out = x;

    function F = eq(X,out2,d2,C2,xdot2,control2,h2)

        % Loop over collocation points
        for jj=1:d2
           % Expression for the state derivative at the collocation point
           xp = C2(1,jj+1)*out2;
           for rr=1:d2
               xp = xp + C2(rr+1,jj+1)*X(:,rr);
           end

           % Append collocation equations
           fj = xdot2(X(:,jj),control2);
           F(:,jj) = h2*fj - xp;

           % Add contribution to the end state
        end
    end
    fun = @(X)eq(X,out,d,C,xdot, u,h);
    pp = xdot(x,u);
    gain = linspace(1,d,d);
    options_fun = optimoptions('fsolve','Display','none');
    x_colloc = fsolve(fun,  x+pp*gain*(h/d), options_fun);% [x+pp*(h/3),x+pp*(2*h/3),x+(pp*h)]);
%     x_next = x_colloc(:,d);
                        
    x_next = D(1)*x;
    x_next = x_next + sum(D(2:end)'.*x_colloc,2);


end