function [best] = run(id,numBest,numGens,options)
arguments
    id
    numBest
    numGens
    options.curve1 = [];
    options.curve2 = [];
    options.stage = 1;
    options.doPlot = true;
end


xcpp1 = XChamberPressureProfile(24,3:5);
xcpp2 = XChamberPressureProfile(24,3);
xcpp1.targetScore = 1;
xcpp2.targetScore = 1;

if ~isempty(options.curve1)
    if numel(options.curve1) == 24 && isa(options.curve1,'double')
        xcpp1 = xcpp1.loadCurve(options.curve1);
    elseif isa(options.curve1, 'XChamberPressureProfile')
        xcpp1 = options.curve1;
    end
end

if ~isempty(options.curve2)
    if numel(options.curve2) == 24 && isa(options.curve2,'double')
        xcpp2 = xcpp2.loadCurve(options.curve2);
    elseif isa(options.curve2, 'XChamberPressureProfile')
        xcpp2 = options.curve2;
    end
end

mf = @(cpp) mutateFun(cpp);

if options.doPlot
    hold off
    figure(1);
    if options.stage == 1
        xcpp1 = mf(xcpp1);
        xcpp1.showPlot();
    elseif option.stage == 2
        xcpp2 = mf(xcpp2);
        xcpp2.showPlot();
    end
    drawnow
    pf = @(best,rest,genNum) plotFun(best,rest,genNum);
else
    pf = @(best,rest,genNum) 0;
end


if options.stage == 1
    ef  = @(cpp) evalFun(cpp,xcpp2,id);
    [best,~] = xcpp1.Trial(ef,mf,numBest,numGens,pf);
elseif options.stage == 2
    ef = @(cpp) evalFun(xcpp1,cpp,id);
    [best,~] = xcpp2.Trial(ef,mf,numGens,pf);
end



ID = best.genNum + 1;
lastnum = 0;
if id > 1
    lastnum = str2double(readmatrix("ARM_Profiles.xls","Range","B2:B127","OutputType","string"));
    lastnum_blankindex = isnan(lastnum);
    lastnum(lastnum_blankindex) = [];
    lastnum = lastnum(end);
end
output_vars = {string(1:length(best.coords))};
output_vars2 = {'First Stage Chamber Pressure Profile'};

out_table = array2table(zeros(1,length(output_vars{1,:})+1),'VariableNames',['Combination',output_vars{1,:}]);
row = 1;
col1 = 1;
col2 = 25;
RangeVariable1 = xlsAddr(row,col1);
RangeVariable2 = xlsAddr(row,col2);
RangeVariable = [RangeVariable1,':',RangeVariable2];
writetable(out_table,'ARM_Profiles.xls','Range',RangeVariable);

import_combo = str2double(readmatrix("Design Matrix - Praeto.csv","Range","A2:A127","OutputType","string"));

output_vec = [import_combo(id),best.coords];
out_table = array2table(output_vec);
row = id + 1;
col1 = 1;
col2 = 25;
RangeVariable1 = xlsAddr(row,col1);
RangeVariable2 = xlsAddr(row,col2);
RangeVariable = [RangeVariable1,':',RangeVariable2];
writematrix(output_vec,'ARM_Profiles.xls','Range',RangeVariable);
end


function out = mutateFun(cpp)
out = cpp + cpp.nudge();
out = out + cpp.yShiftKeyPoints();
out = out + cpp.xShiftKeyPoints();
out = out + cpp.changeKeyPointsNum();
out = out + cpp.reinterpolateCoords();
end


function [score] = evalFun(cpp1,cpp2,iter)
[score] = Atmospheric_1DoF(cpp1,cpp2,iter);
end

function plotFun(best,rest,genNum)
hold off;

best.showPlot(3);
rest.showPlot(1,0.5);
ylim([0,best(1).maxPressure]);
runCount = XChamberPressureProfile.runCounter(best);
title(sprintf("Generation n = %d, Total Runs:%d",genNum,runCount))
legendText = [];
for n = 1:numel(best)
    legendText = [legendText,sprintf("Score = %.6f",best(n).score)];
end

for n = 1:numel(rest)
    legendText = [legendText,sprintf("Score = %.6f",rest(n).score)];
end

legend(legendText)
drawnow();
end