function mass = noseconeMass(shape, shapeParameter, wallThickness, length, AF_thickness, AF_inner_diameter, tipDensity, bodyDensity, tipLength)

    %calc tip volume
    %calc frustum volume
    %multiply by respective material densities
    %add
    
    baseDiameter = 2 * AF_thickness + AF_inner_diameter;
    
    tipVol = 0; % reassigned if there is a solid tip
    baseR = baseDiameter / 2;

    % assign shape equation
    if shape == 1 % power series
        f = @(x) baseR * (x / length) .^ shapeParameter;
    elseif shape == 2 % haak series
        angle = @(x) acos(1 - ((2 * x) / length));
        f = @(x) (baseR / sqrt(pi)) * sqrt(angle(x) - (sin(2 * angle(x)) / 2) + shapeParameter * (sin(angle(x))) .^ 3);
    elseif shape == 3 % conical
        f = @(x) (x * baseR) / length;
    elseif shape == 4 % tangent ogive
        rho = (baseR .^ 2 + length .^ 2) / (2 * baseR);
        f = @(x) sqrt( (rho .^ 2) - (length - x) .^ 2) + baseR - rho;
    elseif shape == 5 %elliptical
        f = @(x) baseR * sqrt(1 - (x .^ 2 / length .^ 2));
    end

    if tipLength > 0
        tipVol = calcTipVol(f, tipLength);
    end

    hollowVol = calcHollowVol(f, wallThickness, tipLength, length);

    tipMass = tipVol * tipDensity;
    hollowMass = hollowVol * bodyDensity;

    mass = tipMass + hollowMass;
end

% CALC VOLUME OF HOLLOW PORTION (FRUSTUM OR CONE)
function hollowVol = calcHollowVol(shapeFunction, wallThickness, tipLength, length)
    % outer f squared - inner f squared inside the integral
    
    outerF = @(x) shapeFunction(x) .^ 2;
    innerF = @(x) (shapeFunction(x) - wallThickness) .^ 2;
    
    integrand = @(x) outerF(x) - innerF(x); 

    hollowVol = pi * integral(integrand, tipLength, length);
end


% CALC VOLUME OF TIP
function tipVol = calcTipVol(shapeFunction, tipLength)
    integrand = @(x) shapeFunction(x) .^2;
    tipVol = pi * integral(integrand, 0, tipLength);
end


