clc
clear
data = csvread('ARM_Metrics.csv');

[index_select,Value_all,Cost_all,comb,mp,mtot,tb1,tb2,maxQ,deltaV,mpl,alt,D1,D2,dvSplit] = Pareto_Selection(data);

% for i = 1:length(alt)
%     figure(1)
%     s = scatter(Cost_all(i),Value_all(i),'bo');
%     row = dataTipTextRow('Combination',comb(i));
%     s.DataTipTemplate.DataTipRows(1) = row;
%     row = dataTipTextRow('Payload Mass',mpl(i));
%     s.DataTipTemplate.DataTipRows(2) = row;
%     row = dataTipTextRow('Propellant Mass',mp(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('Delta V',deltaV(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('Total Mass',mtot(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('1st Stage Burn Time',tb1(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('2nd Stage Burn Time',tb2(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('Max Q',maxQ(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('Altitude',alt(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('1st Stage Diameter',D1(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     row = dataTipTextRow('2nd Stage Diameter',D2(i));
%     s.DataTipTemplate.DataTipRows(end+1) = row;
%     hold on
% end

figure(1)
s = scatter(Cost_all,Value_all,'bo');
row = dataTipTextRow('Combination',comb);
s.DataTipTemplate.DataTipRows(1) = row;
row = dataTipTextRow('Payload Mass',mpl);
s.DataTipTemplate.DataTipRows(2) = row;
row = dataTipTextRow('Propellant Mass',mp);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('Delta V',deltaV);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('Total Mass',mtot);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('1st Stage Burn Time',tb1);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('2nd Stage Burn Time',tb2);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('Max Q',maxQ);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('Altitude',alt);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('1st Stage Diameter',D1);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('2nd Stage Diameter',D2);
s.DataTipTemplate.DataTipRows(end+1) = row;
row = dataTipTextRow('Delta V Split',dvSplit);
s.DataTipTemplate.DataTipRows(end+1) = row;

hold on
s2 = scatter(Cost_all(index_select),Value_all(index_select),'ro');
row2 = dataTipTextRow('Combination',comb(index_select));
s2.DataTipTemplate.DataTipRows(1) = row2;
row2 = dataTipTextRow('Payload Mass',mpl(index_select));
s2.DataTipTemplate.DataTipRows(2) = row2;
row2 = dataTipTextRow('Propellant Mass',mp(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('Delta V',deltaV(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('Total Mass',mtot(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('1st Stage Burn Time',tb1(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('2nd Stage Burn Time',tb2(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('Max Q',maxQ(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('Altitude',alt(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('1st Stage Diameter',D1(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('2nd Stage Diameter',D2(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;
row2 = dataTipTextRow('Delta V Split',dvSplit(index_select));
s2.DataTipTemplate.DataTipRows(end+1) = row2;

xlabel('Cost')
ylabel('Value')
title('Pareto Analysis for Vehicle Sizing')

reddots = data(index_select,:)
