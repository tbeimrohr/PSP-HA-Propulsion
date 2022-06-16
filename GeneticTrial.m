classdef (Abstract) GeneticTrial
    properties
        targetScore = 0
        score
        trialID
        birthGeneration
        scoreFunction
        profile
        allProfiles
%     end
% 
%     properties (Hidden)
        genNum
        trialSetup = false
    end

    methods (Abstract)
        obj = newFrom(obj,numNew);
        iseq = eq(a,b)
    end


    methods


        function obj = set(obj,options)
                propNames = properties(class(obj));
                optionsFields = fields(options);
                for n = 1:numel(optionsFields)
                    if ismember(optionsFields(n),propNames)
                        for m = 1:numel(obj)
                            obj(m).(string(optionsFields(n))) = options.(string(optionsFields(n)));
                        end
                    else
                        error("Unknown property: %s",string(optionsFields(n)));
                    end
                end
        end

%         function obj = set(obj,options)
%             arguments
%                 obj
%                 options.targetScore = nan
%                 options.score = nan
%                 options.birthGeneration = nan
%                 options.scoreFunction = nan
%             end
% 
%             for n = 1:numel(obj)
%                 if ~isnan(options.targetScore)
%                     obj(n).targetScore = options.targetScore;
%                 end
%                 if ~isnan(options.score)
%                     obj(n).score = options.score;
%                 end
%                 if ~isnan(options.birthGeneration)
%                     obj(n).birthGeneration = options.birthGeneration;
%                 end
%                 if ~isnan(options.scoreFunction)
%                     obj(n).scoreFunction = options.scoreFunction;
%                 end
%             end
%         end

        function scores = getScores(obj)
           scores = nan(1,numel(obj));
            for n = 1:numel(obj)
                if isempty(obj(n).score)
                    error("Asked for score not yet calculated.")
                end
                scores(n) = obj(n).score;
            end
        end

        function obj = calculateScore(obj)
            for n = 1:numel(obj)
                if isempty(obj(n).scoreFunction)
                    error("Score function has not been set.")
                end
                if isempty(obj(n).score)
                    GeneticTrial.runCounter(obj(n));
                    obj(n).profile = n;
                    obj(n).score = obj(n).scoreFunction(obj(n));
                end
            end
        end

        function best = getBest(obj,numProfiles,targetScore)
            if isempty(obj(1).scoreFunction)
                error("Score Function has not been specified");
            end

            if nargin < 3
                if isempty(obj(1).targetScore) 
                    error("Target Score has not been specified");
                else
                    targetScore = obj.targetScore;
                end

            end
            if nargin < 2
                numProfiles = 1;
            end

            scores = obj.getScores();

            if targetScore == inf
                targetScore = max(scores);
            elseif targetScore == -inf
                targetScore = min(scores);
            end

            sortDirection = 'ascend';
            if numProfiles < 0
                numProfiles = abs(numProfiles);
                sortDirection = 'descend';
            end

            [~,sortIdx] = sort(abs(targetScore-scores),sortDirection);
            sortedScores = scores(sortIdx);
            
%             if strcmp(sortDirection,'ascend')
                best = obj(sortIdx(1:numProfiles));
%                 best = obj(scores <= sortedScores(numProfiles));
%             else
%                 best = obj(scores >= sortedScores(numProfiles));
%             end
            if numel(best) > numProfiles
                best = best(1:numProfiles);
            end
%             n = numProfiles
%             
        end

        function obj = setBirthGeneration(obj,value)
            for n = 1:numel(obj)
                if isempty(obj(n).birthGeneration)
                    obj(n).birthGeneration = value;
                end
                obj(n).genNum = value;
            end
        end

        function obj = plus(obj1, obj2)
            obj = [obj1 obj2];
        end

        function obj = setupTrial(obj,evalFun)
            persistent trialID
            if isempty(trialID)
                trialID = 0;
            else
                trialID = trialID + 1;
            end

            options.trialID = trialID;
            options.birthGeneration = 0;
            options.scoreFunction = evalFun;
            options.trialSetup = true;
            options.genNum = 0;
            options.profile = [];
            options.allProfiles = numel(obj);
            obj = obj.set(options);
            
            val = GeneticTrial.runCounter(obj,0);

        end

        function [best,rest,newGen] = TrialGeneration(obj,mutateFun,numBest,plotFun)
            if nargin < 4
                plotFun = @(best,rest,genNum) 0;
            end
            
            for n = 1:numel(obj)
                if obj(n).trialSetup == false
                    error("Trial not setup.");
                end
%                 obj(n).info.num = n;
%                 obj(n).info.total = numel(obj);
            end

            obj = obj.calculateScore();
%             numel(obj)
%             numBest
            rest = obj.getBest(numBest-numel(obj));
            best = obj.getBest(numBest);
            
            plotFun(best,rest,best(1).genNum);

            newGen = mutateFun(best);
            newGen = newGen.setBirthGeneration(obj(1).genNum+1);

        end

        
        function [best,rest,ID] = Trial(obj,evalFun,mutateFun,numBest,numGenerations,plotFun)
            if nargin < 6
                plotFun = @(best,rest,genNum) 0;
            end
            
           newGen = obj.setupTrial(evalFun);

            for n = 0:numGenerations
                if n > 0
                    import_metrics = readmatrix("ARM_Metrics_temp.xls");
                    select_metrics = import_metrics(find(import_metrics(:,17) == max(import_metrics(:,17))),:);
                    row = 2;
                    col1 = 1;
                    col2 = 18;
                    RangeVariable1 = xlsAddr(row,col1);
                    RangeVariable2 = xlsAddr(row,col2);
                    RangeVariable = [RangeVariable1,':',RangeVariable2];
                    writematrix(select_metrics,'ARM_Metrics_temp.xls','Range',RangeVariable);

                    import_profiles = readmatrix("ARM_Profiles_temp.xls");
                    select_profile = import_profiles(find(import_metrics(:,17) == max(import_metrics(:,17))),:);
                    row = 1;
                    col1 = 1;
                    col2 = 24;
                    RangeVariable1 = xlsAddr(row,col1);
                    RangeVariable2 = xlsAddr(row,col2);
                    RangeVariable = [RangeVariable1,':',RangeVariable2];
                    writematrix(select_profile,'ARM_Profiles_temp.xls','Range',RangeVariable)
                end
                [best,rest,newGen] = newGen.TrialGeneration(mutateFun,numBest,plotFun);
            end

            ID = obj.trialID;

        end

    end

    
    methods (Static)    
    
     function runNumber = runCounter(obj,setValue)
            persistent counter
%             key = sprintf("%s:%s",class(obj),func2str(obj.scoreFunction));
            if numel(obj) > 1
                obj = obj(1);
            end
            
            key = string(obj.trialID);
            if isempty(counter)
                counter = containers.Map;
            end
            if nargin == 2
                counter(key) = setValue;
            elseif nargout == 0
                counter(key) = counter(key) + 1;
            end
            runNumber = counter(key);

        end
    end
end
