function     [delH1,delH2] = computeDSLHeatUnary(P, y1, T, PRef, TRef, parameters)
R = 8.3145;
% Analytical computation of heat of adsorption as a function of P and T
    b1 = parameters.bo_1 .* exp(-parameters.delUb_1 ./ (R .* T.*TRef)); % Equilibrium constant for 1 on site b [m3/mol]
    d1 = parameters.do_1 .* exp(-parameters.delUd_1 ./ (R .* T.*TRef)); % Equilibrium constant for 1 on site d [m3/mol]
    b2 = parameters.bo_2 .* exp(-parameters.delUb_2 ./ (R .* T.*TRef)); % Equilibrium constant for 2 on site b [m3/mol]
    d2 = parameters.do_2 .* exp(-parameters.delUd_2 ./ (R .* T.*TRef)); % Equilibrium constant for 2 on site d [m3/mol]

    c1 =     y1.*P.*PRef./(R.*T.*TRef); % Bulk concentration of 1 [mol/m3]
    c2 = (1-y1).*P.*PRef./(R.*T.*TRef); % Bulk concentration of 2 [mol/m3]

    delH1 = -(parameters.qsb_1.*b1.*parameters.delUb_1./(1+b1.*c1).^2 + parameters.qsd_1.*d1.*parameters.delUd_1./(1+d1.*c1).^2)./...
        (parameters.qsb_1.*b1./((1+b1.*c1).^2)+parameters.qsd_1.*d1./((1+d1.*c1).^2)); % RT^2 * dlnp/dT (Clausius Clapeyron)

    delH2 = -(parameters.qsb_2.*b2.*parameters.delUb_2./((1+b2.*c2).^2) + parameters.qsd_2.*d2.*parameters.delUd_2./((1+d2.*c2).^2))./...
        (parameters.qsb_2.*b2./((1+b2.*c2).^2)+parameters.qsd_2.*d2./((1+d2.*c2).^2)); % RT^2 * dlnp/dT (Clausius Clapeyron)
end