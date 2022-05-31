function M = NewtonsMethod_P1(mdot,mdot_tot,p0,Ap,R,Tc,gamma,g,M_guess)
%% Newton's Method Solver Eq. 4.24
% This function solves equation 4.24 in the SRM section for Mach

% All inputs should be in imperial units; If in SI, the equation will
% change (there will not be a g in the sqrt)
% The guess for mach should be in the range of 0 - 0.5 
% (I'm guessing mach from previous axial step as the seed)

error = 1;
while error > 0.001 
    M_save = M_guess;
    
    a = g*p0*Ap / (3.281*(sqrt(R*Tc/(gamma)))); %front coefficient 3.281
    b = (gamma - 1) / 2; %inner coefficient
    c = -(gamma+1)/(2*(gamma-1)); %exponent
    
    fm = mdot_tot+mdot/2 - a*M_guess*(1+b*M_guess^2)^c; %original function
    fm_prime = -a * ((1+b*M_guess^2)^c + 2*b*c*M_guess^2*(1+b*M_guess^2)^(c-1)); %derivative
        
    M_guess = M_guess - fm/fm_prime; %updating guess
    error = abs(M_guess - M_save); %comparing new and old guesses to see if they meet convergence criterion
end
    
M = M_guess;
