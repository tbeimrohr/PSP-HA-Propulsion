function [index_select,Value_all,Cost_all,comb,mp,mtot,tb1,tb2,maxQ,deltaV,mpl,alt,D1,D2,dvSplit] = Pareto_Selection(data)

%% Gather Data
comb = data(:,1);
mp = data(:,2);
mi = data(:,3);
mtot = data(:,4);
tb1 = data(:,5);
tb2 = data(:,6);
maxQ = data(:,7);
deltaV = data(:,8);
mpl = data(:,9);
alt = data(:,10);
D1 = data(:,11);
D2 = data(:,12);
dvSplit = data(:,13);

for i = 1:length(comb)
    [tb_factor(i),mpl_factor(i),alt_factor(i),mp_factor(i),maxQ_factor(i),L_D_factor(i),L_D(i)] = Pareto(mpl(i),D1(i),D2(i),alt(i),mp(i),tb1(i),tb2(i),maxQ(i));
end

%% Normalize
power = 1;
tb_factor = (tb_factor./max(tb_factor));
mpl_factor = (mpl_factor./max(mpl_factor)).^power;
alt_factor = (alt_factor/max(alt_factor));
mp_factor = (mp_factor./max(mp_factor));
maxQ_factor = (maxQ_factor./max(maxQ_factor));
L_D_factor = (L_D_factor./max(L_D_factor));


%% Weights

% Weight 1
f1 = .4; %tb weight                     bad
f2 = .05; % payload weight              good
f3 = .35; %altitude weight               good
f4 = .4; % propellant mass weight       bad
f5 = .2; % max Q weight                 bad
f6 = .6; %L/D weight                   good

% % Weight 2
% f1 = .1;
% f2 = .05;
% f3 = .65;
% f4 = .6;
% f5 = .3;
% f6 = .3;

% Weight Group
% f1 = 0.125;
% f2 = 0.039;
% f3 = 0.543;
% f4 = 0.479;
% f5 = 0.396;
% f6 = 0.418;

%% Value and Cost

Value_all = mpl_factor.*f2 + alt_factor.*f3 + L_D_factor.*f6;
Cost_all = tb_factor.*f1 + mp_factor.*f4 + maxQ_factor.*f5;

desired = [0 (max(Value_all)+mean(Value_all))/1.75];
for i = 1:length(comb)
    dist(i) = distance(desired,[Cost_all(i) Value_all(i)]);
end
selection1 = ((dist))./max(dist);
min_value = min(selection1);
mean_value = mean(selection1);
index_select = find(selection1 >= min_value & selection1 <= .9*mean_value);

