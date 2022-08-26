function [tb_factor,mpl_factor,alt_factor,mp_factor,maxQ_factor,L_D_factor,L_D] = Pareto(mpl,D1,D2,alt,mp,tb1,tb2,maxQ)

%% Referance Values (USC)
tb_ref = 14; %sec
alt_ref = 103.57; %km
L_D_ref = 19.5;
mass_tot = 140.614; %kg
maxQ_ref = 570; %kpa
mpl_ref = 10; %kg

%% Constraints
rho = 0.0686165657; % lbm/in3
V_grain = 32.12955906; %in3
V_loading = .8;



%% Calculations
tb_factor = (tb1+tb2)/(tb_ref);
mpl_factor = mpl/mpl_ref;
alt_factor = alt/alt_ref;
mp_factor = mp/mass_tot;
maxQ_factor = maxQ/maxQ_ref;

L_D = ((mp*2.205/V_loading)/rho)/(pi*(((D1+D2)/4)^2))/((D1+D2)/2);
L_D_factor = 1 - ((L_D_ref - L_D)/L_D)^2;

