function [index_select,Value_all,Cost_all,comb,mp,mtot,tb1,tb2,maxQ,deltaV,mpl,alt,D1,D2] = Pareto_Selection(data)

%% Gather Data
try
    comb = data{:,1};
    mp = data{:,2};
    mi = data{:,3};
    mtot = data{:,4};
    tb1 = data{:,5};
    tb2 = data{:,6};
    maxQ = data{:,7};
    deltaV = data{:,8};
    mpl = data{:,9};
    alt = data{:,10};
    D1 = data{:,11};
    D2 = data{:,12};
catch
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
end

for i = 1:length(comb)
    [tb_factor(i),mpl_factor(i),alt_factor(i),mp_factor(i),maxQ_factor(i),L_D_factor(i)] = Pareto(mpl(i),D1(i),D2(i),alt(i),mp(i),tb1(i),tb2(i),maxQ(i));
end

%% Normalize
% P = 3;
% tb_factor = (tb_factor./max(tb_factor)).^P;
% mpl_factor = (mpl_factor./max(mpl_factor)).^P;
% alt_factor = (alt_factor/max(alt_factor)).^P;
% mp_factor = (mp_factor./max(mp_factor)).^P;
% maxQ_factor = (maxQ_factor./max(maxQ_factor)).^P;
% % L_D_factor = (L_D_factor ./ max(L_D_factor)).^4;

%% Weights

% Weight 1
% f1 = .2; % time - bad
% f2 = .05; % payload - good
% f3 = .5; % altitude - good
% f4 = .4; % propellant mass factor - bad
% f5 = .4; % dynamic pressure - bad
% f6 = .45; % LD Factor - good

% % Weight 2
% f1 = .125;
% f2 = .05;
% f3 = .65;
% f4 = .6;
% f5 = .275;
% f6 = .3;

% Weight 3
f1 = .2;
f2 = .1;
f3 = .4;
f4 = .1;
f5 = .7;
f6 = .5;

%% Value and Cost

Value_all = mpl_factor.*f2 + alt_factor.*f3 + L_D_factor.*f6;
Cost_all = tb_factor.*f1 + mp_factor.*f4 + maxQ_factor.*f5;

desired = [0 max(Value_all)];
for i = 1:length(comb)
    dist(i) = distance(desired,[Cost_all(i) Value_all(i)]);
end
selection1 = ((dist))./max(dist);
min_value = min(selection1);
mean_value = mean(selection1);
index_select = find(selection1 >= min_value & selection1 <= 1.025*mean_value);
