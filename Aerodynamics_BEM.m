% -----------------------------
% Function: Calculates aerodynamic torque and thrust based on the "Blade
% Element Momentum Theory following [3]
% ------------
% Input:
% - x           vector of states 
% - v_0         scalar of rotor-effective wind speed
% - Parameter   struct of Parameters
% ------------
% Output:
% - M_a         scalar of aerodynamic torque
% - F_a         scalar of aerodynamic thrust
% ----------------------------------
function [M_a,F_a] = Aerodynamics_BEM(x,v_0,Parameter)

BEM  = Parameter.BEM;
rho  = Parameter.General.rho;
Omega = x(1);

M_a = 0;
F_a = 0;

B = 3; % number of blades (NREL 5MW)
R_tip = 63; % total rotor radius

for i = 1:BEM.BldNodes

    r   = BEM.RNodes(i);
    c   = BEM.Chord(i);
    dr  = BEM.DRNodes(i);
    foil = BEM.NFoil(i);

    % initial guesses
    a  = 0.3;
    ap = 0.0;

    for iter = 1:50

        % local velocities
        Vax  = v_0 * (1 - a);
        Vtan = Omega * r * (1 + ap);

        Vrel = sqrt(Vax^2 + Vtan^2);
        phi  = atan2(Vax, Vtan);

        % angle of attack (deg!)
        alpha = phi * 180/pi - BEM.AeroTwst(i);

        % airfoil data
        CL = interp1(BEM.AoA{foil}, BEM.Cl{foil}, alpha, 'linear', 'extrap');
        CD = interp1(BEM.AoA{foil}, BEM.Cd{foil}, alpha, 'linear', 'extrap');

        % force coefficients (Cn and Ct modified to include drag)
        Cn = CL*cos(phi) + CD*sin(phi);
        Ct = CL*sin(phi) - CD*cos(phi);

        % Prandtl tip-loss factor calculation
        if sin(phi) > 0.001
            f_tip = (B / 2) * ((R_tip - r) / (r * sin(phi)));
            F = (2 / pi) * acos(exp(-f_tip));
            F = max(F, 1e-4); % Avoid division by zero
        else
            F = 1;
        end

        % solidity
        sigma = (B * c) / (2*pi*r);

        % induction updates (including Prandtl loss factor F)
        a_new  = 1 / ((4*F*sin(phi)^2)/(sigma*Cn) + 1);
        ap_new = 1 / ((4*F*sin(phi)*cos(phi))/(sigma*Ct) - 1);

        % relaxation
        a  = 0.75*a + 0.25*a_new;
        ap = 0.75*ap + 0.25*ap_new;

        % convergence check
        if abs(a-a_new) < 1e-4 && abs(ap-ap_new) < 1e-4
            break;
        end
    end

    % =========================
    %  BEM FORCE MODEL
    % =========================

    dT = 0.5 * rho * Vrel^2 * Cn * c * dr;
    dQ = 0.5 * rho * Vrel^2 * Ct * c * r * dr;

    % sum over blades
    F_a = F_a + B * dT;
    M_a = M_a + B * dQ;

end

end
