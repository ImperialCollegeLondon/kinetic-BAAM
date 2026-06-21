function     [delH1,delH2] = computeQHeatUnary(P, y1, T, q1, q2, PRef, TRef, qRef, parameters)
% R = 8.3145;
% Analytical computation of heat of adsorption as a function of P and T
    % b1 = parameters.bo_1 .* exp(-parameters.delUb_1 ./ (R .* T.*TRef)); % Equilibrium constant for 1 on site b [m3/mol]
    % d1 = parameters.do_1 .* exp(-parameters.delUd_1 ./ (R .* T.*TRef)); % Equilibrium constant for 1 on site d [m3/mol]
    % b2 = parameters.bo_2 .* exp(-parameters.delUb_2 ./ (R .* T.*TRef)); % Equilibrium constant for 2 on site b [m3/mol]
    % d2 = parameters.do_2 .* exp(-parameters.delUd_2 ./ (R .* T.*TRef)); % Equilibrium constant for 2 on site d [m3/mol]

%     c10 =     y1.*P.*PRef./(R.*T.*TRef); % Bulk concentration of 1 [mol/m3]
%     c20 = (1-y1).*P.*PRef./(R.*T.*TRef); % Bulk concentration of 2 [mol/m3]
% 
%     q1_fun = @(c1) sqrt((q1.*qRef-((parameters.qsb_1 .* b1 .* c1 ) ./   (1 + b1 .* c1 ) + parameters.qsd_1 .* d1 .* c1 ./ (1 + (d1  .* c1) ))).^2) ;%mol/kg
%     q2_fun = @(c2) sqrt((q2.*qRef-((parameters.qsb_2 .* b2 .* c2 ) ./   (1 + b2 .* c2 ) + parameters.qsd_2 .* d2 .* c2 ./ (1 + (d2  .* c2) ))).^2) ;%mol/kg
%     options = optimoptions('fmincon', 'Display', 'off');
%     c1 = fmincon(q1_fun,c10,[],[],[],[],[],[],[],options);
%     c2 = fmincon(q2_fun,c20,[],[],[],[],[],[],[],options);
% 
%     delH1 = -(parameters.qsb_1.*b1.*parameters.delUb_1./(1+b1.*c1).^2 + parameters.qsd_1.*d1.*parameters.delUd_1./(1+d1.*c1).^2)./...
%         (parameters.qsb_1.*b1./((1+b1.*c1).^2)+parameters.qsd_1.*d1./((1+d1.*c1).^2)); % RT^2 * dlnp/dT (Clausius Clapeyron)
% 
%     delH2 = -(parameters.qsb_2.*b2.*parameters.delUb_2./((1+b2.*c2).^2) + parameters.qsd_2.*d2.*parameters.delUd_2./((1+d2.*c2).^2))./...
%         (parameters.qsb_2.*b2./((1+b2.*c2).^2)+parameters.qsd_2.*d2./((1+d2.*c2).^2)); % RT^2 * dlnp/dT (Clausius Clapeyron)

      q1 = max(q1,eps);
      q2 = max(q2,eps);
      delH1 = interp1(parameters.heat1(:,1),parameters.heat1(:,2),q1.*qRef);  
      delH2 = interp1(parameters.heat2(:,1),parameters.heat2(:,2),q2.*qRef);  


end