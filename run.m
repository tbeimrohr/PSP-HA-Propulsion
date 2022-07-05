function [best] = run(id,numBest,numGens,options)
arguments
    id
    numBest
    numGens
    options.curve1 = [];
    options.curve2 = [];
    options.stage = [];
    options.doPlot = true;
end


xcpp1 = XChamberPressureProfile(24,[3 4]);
xcpp2 = XChamberPressureProfile(24,[3 4]);
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
    if options.stage == 1
        figure(1);
        xcpp1 = mf(xcpp1);
        xcpp1.showPlot();
        xlabel('Index')
        ylabel('Pressure [psi]')
        title(sprintf("Chamber Pressure Profile"))
        txt = sprintf('img1.png');
        saveas(1,txt)
    elseif options.stage == 2
        figure(1);
        xcpp2 = mf(xcpp2);
        xcpp2.showPlot();
        xlabel('Index')
        ylabel('Pressure [psi]')
        title(sprintf("Chamber Pressure Profile"))

        txt = sprintf('img2.png');
        saveas(1,txt)
    elseif options.stage == [1 2]
        figure(1);
        xcpp1 = mf(xcpp1);
        xcpp1.showPlot();
        xlabel('Index')
        ylabel('Pressure [psi]')
        figure(2);
        xcpp2 = mf(xcpp2);
        xcpp2.showPlot();
        xlabel('Index')
        ylabel('Pressure [psi]')
    end

    drawnow
    pf = @(best,rest,genNum) plotFun(best,rest,genNum);
else
    pf = @(best,rest,genNum) 0;
end


if options.stage == 1
    ef  = @(cpp) evalFun(cpp,xcpp2,id);
    [best] = xcpp1.Trial(ef,mf,numBest,numGens,pf);
elseif options.stage == 2
    ef = @(cpp) evalFun(xcpp1,cpp,id);
    [best,~] = xcpp2.Trial(ef,mf,numBest,numGens,pf);
elseif options.stage == [1 2]
    ef  = @(cpp) evalFun(cpp,xcpp2,id);
    [best,~] = xcpp1.Trial(ef,mf,numBest,numGens,pf);

    ef = @(cpp) evalFun(xcpp1,cpp,id);
    [best2,~] = xcpp2.Trial(ef,mf,numBest,numGens,pf);
end


if length(options.stage) == 1
    %% metrics

    choice = best(1,1);
    import_combo = fileread("combos.txt");
    temporary = convertCharsToStrings(import_combo);
    iter = (strfind(temporary,id)-1)/14 + 1;

    output_vars = {'Combination','Propellant Mass [kg]','Inert Mass [kg]',...
        'Total Mass [kg] (including payload)','First Stage Burn Time [sec]',...
        'Second Stage Burn Time [sec]','Maximum Dynamic Pressure (Max Q) [kPa]',...
        'Delta V [km/s]','Payload Mass (kg)','Alt [km]','Diameter 1','Diameter 2','Delta V Split for 1st Stage','Coast limit','Value','Cost','Score'};


    out_table = array2table(zeros(1,length(output_vars)), 'VariableNames',output_vars);
    row = 1;
    col1 = 1;
    col2 = 17;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    writetable(out_table,'ARM_Metrics.xls','Range',RangeVariable);

    output_vec = readmatrix("ARM_Metrics_temp.xls");
    choice_profile = find(output_vec(:,17) == best(1,1).score);
    choice_metrics = output_vec(choice_profile,1:17);

    row = iter+1;
    col1 = 1;
    col2 = 17;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    writematrix(choice_metrics,'ARM_Metrics.xls','Range',RangeVariable);

    delete('ARM_Metrics_temp.xls')

    %% profiles

    %     ID = best.genNum + 1;

    if length(best) > 1
        for n = 1:numBest
            score(n) = best(1,n).score;
        end
        select = find(score == max(score));
        choice = best(1,select);
    end
    %     lastnum = 0;
    %     if id > 1
    %         lastnum = str2double(readmatrix("ARM_Profiles.xls","Range","B2:B127","OutputType","string"));
    %         lastnum_blankindex = isnan(lastnum);
    %         lastnum(lastnum_blankindex) = [];
    %         lastnum = lastnum(end);
    %     end
    output_vars = {string(1:length(choice.coords))};
    %     output_vars2 = {'First Stage Chamber Pressure Profile'};

    out_table = array2table(zeros(1,length(output_vars{1,:})+1),'VariableNames',['Combination',output_vars{1,:}]);
    row = 1;
    col1 = 1;
    col2 = 25;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    if options.stage == 1
        name = "ARM_Profiles.xls";
    elseif options.stage == 2
        name = "ARM_Profiles2.xls";
    end
    writetable(out_table,name,'Range',RangeVariable);

    %     import_combo = str2double(readmatrix("Design Matrix - Pareto.csv","Range","A2:A127","OutputType","string"));

    output_vec = [str2num(id),choice.coords];
    out_table = array2table(output_vec);

    row = iter + 1;
    col1 = 1;
    col2 = 25;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    writematrix(output_vec,name,'Range',RangeVariable);
    delete('ARM_Profiles_temp.xls')
else
    ID = best.genNum + 1;
    lastnum = 0;
    %     if id > 1
    %         lastnum = str2double(readmatrix("ARM_Profiles.xls","Range","B2:B127","OutputType","string"));
    %         lastnum_blankindex = isnan(lastnum);
    %         lastnum(lastnum_blankindex) = [];
    %         lastnum = lastnum(end);
    %     end
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

    %     import_combo = str2double(readmatrix("Design Matrix - Praeto.csv","Range","A2:A127","OutputType","string"));

    output_vec = [import_combo(id),best.coords];
    out_table = array2table(output_vec);
    row = id + 1;
    col1 = 1;
    col2 = 25;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    writematrix(output_vec,'ARM_Profiles.xls','Range',RangeVariable);

    %%
    ID = best2.genNum + 1;
    lastnum = 0;
    %     if id > 1
    %         lastnum = str2double(readmatrix("ARM_Profiles2.xls","Range","B2:B127","OutputType","string"));
    %         lastnum_blankindex = isnan(lastnum);
    %         lastnum(lastnum_blankindex) = [];
    %         lastnum = lastnum(end);
    %     end
    output_vars = {string(1:length(best2.coords))};
    output_vars2 = {'Second Stage Chamber Pressure Profile'};

    out_table = array2table(zeros(1,length(output_vars{1,:})+1),'VariableNames',['Combination',output_vars{1,:}]);
    row = 1;
    col1 = 1;
    col2 = 25;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    writetable(out_table,'ARM_Profiles2.xls','Range',RangeVariable);

    import_combo = str2double(readmatrix("Design Matrix - Praeto.csv","Range","A2:A127","OutputType","string"));

    output_vec = [import_combo(id),best2.coords];
    out_table = array2table(output_vec);
    row = id + 1;
    col1 = 1;
    col2 = 25;
    RangeVariable1 = xlsAddr(row,col1);
    RangeVariable2 = xlsAddr(row,col2);
    RangeVariable = [RangeVariable1,':',RangeVariable2];
    writematrix(output_vec,'ARM_Profiles.xls','Range',RangeVariable);
end
end


function out = mutateFun(cpp)
placeholder = cpp;
placeholder.coords = smooth(cpp.coords)';
out = cpp + placeholder;
out = out + cpp.nudge().yShiftKeyPoints().xShiftKeyPoints().changeKeyPointsNum().reinterpolateCoords().reinterpolateCoords();
out = out + cpp.yShiftKeyPoints().changeKeyPointsNum().reinterpolateCoords().reinterpolateCoords().reinterpolateCoords();
out = out + cpp.nudge().nudge().xShiftKeyPoints().changeKeyPointsNum().reinterpolateCoords();
out = out + cpp.nudge().changeKeyPointsNum().reinterpolateCoords();
out = out + cpp.nudge().yShiftKeyPoints().yShiftKeyPoints().xShiftKeyPoints().changeKeyPointsNum().reinterpolateCoords();
out = out + cpp.nudge().xShiftKeyPoints().changeKeyPointsNum();
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
title(sprintf("Chamber Pressure Profile\nGeneration n = %d, Total Runs:%d",genNum,runCount))
xlabel('Index')
ylabel('Pressure [psi]')
legendText = [];
for n = 1:numel(best)
    legendText = [legendText,sprintf("Score = %.3f",best(n).score)];
end

for n = 1:numel(rest)
    legendText = [legendText,sprintf("Score = %.3f",rest(n).score)];
end

legend(legendText,'location','best')
txt = sprintf('img%d.png',best.trialID+best.genNum);
saveas(1,txt)
drawnow();
end