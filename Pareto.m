function [tb_factor,mpl_factor,alt_factor,mp_factor,maxQ_factor,L_D_factor] = Pareto(mpl,D1,D2,alt,mp,tb1,tb2,maxQ)

%% Referance Values (USC)
tb_ref = 14; %sec
alt_ref = 103.57; %km
L_D_ref = 19.5;
mass_tot = 140.614; %kg
mp_ref = mass_tot*.85;
maxQ_ref = 288; %kpa
mpl_ref = 10; %kg

%% Constraints
rho = 0.0686165657; % lbm/in3
V_grain = 32.12955906; %in3
V_loading = .8;



%% Calculations
tb = tb1 + tb2;
% tb_factor = (tb1+tb2)/(tb_ref);
% mpl_factor = mpl/mpl_ref;
% alt_factor = alt/alt_ref;
% mp_factor = mp/mass_tot;
% maxQ_factor = maxQ/maxQ_ref;
% 
L_D = ((mp*2.205/V_loading)/rho)/(pi*(((D1+D2)/4)^2))/((D1+D2)/2);
L_D_factor = 1 - (abs((L_D_ref - L_D)/L_D))^2;

% 
% tb_factor = (abs((tb_ref - (tb1+tb2))/(tb1+tb2)))/tb_ref;
% mpl_factor = (abs((mpl_ref - (mpl))/(mpl)))/mpl_ref;
% alt_factor = (abs((alt_ref - (alt))/(alt)))/alt_ref;
% mp_factor = (abs((mass_tot - (mp))/(mp)))/mass_tot;
% maxQ_factor = (abs((maxQ_ref - (maxQ))/(maxQ)))/maxQ_ref;
% 
% L_D_factor = 1 - ((abs((L_D_ref - L_D)/L_D))/L_D_ref)^.25
% 
tb_factor = 1 - tb_ref/tb;
mpl_factor = mpl/(mpl_ref);
alt_factor = alt/alt_ref - 1;
mp_factor = 1 - (mp_ref - (mp))/(mp_ref);
maxQ_factor = 1 - maxQ_ref/maxQ;

if tb_factor > 1 
    tb_factor = 1;
elseif tb_factor < 0
    tb_factor = 0;
end
if mpl_factor < 0
    mpl_factor = 0;
end
if alt_factor < 0
    alt_factor = 0;
elseif alt_factor > 1
    alt_factor = 1;
end
if mp_factor < 0
    mp_factor = 0;
end
if maxQ_factor < 0
    maxQ_factor = 0;
end





