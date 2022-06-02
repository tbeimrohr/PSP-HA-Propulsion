classdef XChamberPressureProfile < GeneticTrial
% XChamberPressureProfile Class for genetic optimization of curves
%   Class for generating pressure curves and refining them using trials
%   using a genetic algorithm. Curves are interpolated between a given
%   number of "Key Points"
    properties
        maxPressure = 1750          % maximum chamber pressure
        minStartPercentage = 0.2    % minimum starting % of maxPressure 
        maxStartPercentage = 0.8    % maximum starting % of maxPressure    
        maxPeakRatio = 0.85;         % max height of subsequent peaks
        xShiftConstant = 0.1;       % +/- allowed shift to next key point in xShift
        yShiftConstant = 0.33;       % +/- percent range for yShift mutation
        yIsInt = true;              % forces pressure to be rounded ints (is faster)
        showKeyPoints = false       % displays key points when showPlot() is called
        
        coords              % array storing pressures
        keyPoints           % matrix storing key point values and locations
        numKeyPointsDomain  % array containing allowed number of key points
    end

    methods
        function obj = XChamberPressureProfile(numCoords,numKeyPoints,numProfiles)
            if nargin < 3
                numProfiles = 1;
            end
            obj.numKeyPointsDomain = numKeyPoints;
            if numel(numKeyPoints) > 1
                numKeyPoints = numKeyPoints(randi([1,numel(numKeyPoints)]));
            end
            
            % Evenly distribute key point x values
            obj.keyPoints(1,:) = round(linspace(1,numCoords,numKeyPoints+2));

            % Set first key point within range of start percentages
            obj.keyPoints(2,1) = obj.randInRange(obj.minStartPercentage*obj.maxPressure,obj.maxStartPercentage*obj.maxPressure);
            
            % Set second key point above first key point and below max
            obj.keyPoints(2,2) = obj.randInRange(obj.keyPoints(2,1),obj.maxPressure);
            
            % Set remaining key points with a lower bound of zero, and an
            % upper bound of either double the previous value or the
            % maximum allowed value as calculated by maxNextPeak()
            for n = 3:numKeyPoints+1
                obj.keyPoints(2,n) = obj.randInRange(0,min(obj.maxNextPeak(n),obj.keyPoints(2,n-1)*2));
            end
            
            % Initialize coords and interpolate between key points
            obj.coords = zeros(1,obj.keyPoints(1,end));
            obj = obj.interpolateCoords();
            
            % Repeat initialization for specified number of profiles
            if numProfiles > 1
                obj = obj + XChamberPressureProfile(numCoords,obj.numKeyPointsDomain,numProfiles-1);
            end
        end

        function newProfiles = newFrom(obj,numProfiles)
        % newFrom Creates random new profile(s)
        %   Creates new profile(s) with no resemblance to the parent other
        %   than the same number of coords and numKeyPointsDomain.
        %
        % newProfiles = obj.newFrom() creates a single new profile.
        %
        % newProfiles = obj.newFrom(numProfiles) creates the specified
        %   number of profiles.
        
            if numel(obj) > 1
                obj = obj(1);
            end
            if nargin < 2
                numProfiles = 1;
            end
            newProfiles = XChamberPressureProfile(numel(obj.coords),obj.numKeyPointsDomain,numProfiles);
            
            options.targetScore = obj.targetScore;
            options.trialID = obj.trialID;
            options.scoreFunction = obj.scoreFunction;
            options.trialSetup = obj.trialSetup;
            newProfiles = newProfiles.set(options);
        end

        function obj = loadCurve(obj,curve)
            obj.coords = curve;
            for n = 1:size(obj.keyPoints,2)
                obj.keyPoints(2,n) = curve(obj.keyPoints(1,n));
            end
        end


        function maxPeak = maxNextPeak(obj,nKP)
        % maxNextPeak determines maximum allowed next peak
        %   Determines maximum allowed next peak as a function of the
        %   previous peak height and the specified maxPeakRatio. If there
        %   are no previous peak, the max allowed value is maxPressure.
        %
        % maxPeak = obj.maxNextPeak(nKP) returns the max allowed peak for
        %   the key point with index nKP.

            negSeg = @(n) obj.keyPoints(2,n) < obj.keyPoints(2,n-1);
            if nKP <= 3
                maxPeak = obj.maxPressure;
            elseif negSeg(nKP-1) && ~negSeg(nKP-2)
                maxPeak = obj.keyPoints(2,nKP-2) * obj.maxPeakRatio;
            else
                maxPeak = obj.maxNextPeak(nKP-1);
            end
        end
        
        function obj = xShiftKeyPoints(obj, numProfiles)
        % xShiftKeyPoints randomly redistributes key points horizontally
        %   Acts as a function wrapper to obj.single_xShiftKeyPoints
        %   Function wrapper is required to handle calls when the parent
        %   object is an array (multiplicity > 1) or when multiple profiles
        %   need to be generated.
        %
        % obj = obj.xShiftKeyPoints() creates a "x shifted" profile for 
        %   each element in obj.
        %
        % obj = obj.xShiftKeyPoints(numProfiles) creates the specified
        %   number of "x shifted" profiles for each element in obj.
        %
        %   See also single_xShiftKeyPoints
            if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n) = obj(n).single_xShiftKeyPoints();
            end

            if numProfiles > 1
                obj = obj + obj.xShiftKeyPoints(numProfiles-1);
            end
        end

        function obj = single_xShiftKeyPoints(obj)
        % single_xShiftKeyPoints implements x shifting key points 
        %   The "single" version of this function is the actual x shift
        %   implementation for when object multiplicity is 1. Higher 
        %   multiplicity calls are handled by the above wrapper function.
        %
        % See also xShiftKeyPoints

            r = obj.xShiftConstant;

            % Anonymous function calculates bounds each key point can shift
            calcCoordBound = @(n,r) round(obj.keyPoints(1,n)*r+obj.keyPoints(1,n+1)*(1-r));
            
            % Calculate shift for all key points (excluding first and last)
            for n = 2:size(obj.keyPoints,2)-1
                % Calculate x shift amount
                newX = randi([calcCoordBound(n-1,r),calcCoordBound(n,1-r)]);
                dx = newX - obj.keyPoints(1,n);
                
                % Anonymous function generates coord addresses for range
                coord = @(n)  obj.coords(obj.keyPoints(1,n-1):obj.keyPoints(1,n));

                % Resize coord ranges around new key point location
                leftCoords = imresize(coord(n),[1,numel(coord(n))+dx],"bilinear");
                rightCoords = imresize(coord(n+1),[1,numel(coord(n+1))-dx],"bilinear");
                
                % Update coords and key points with new values
                newCoords = [leftCoords(1:end-1),rightCoords];
                obj.coords(obj.keyPoints(1,n-1):obj.keyPoints(1,n+1)) = newCoords;
                obj.keyPoints(1,n) = newX;
                obj.coords(newX) = obj.keyPoints(2,n);
            end

            % Fix weird artifacts idr exactly
            for n = 1:2
                obj.coords(obj.keyPoints(1,n)) = obj.keyPoints(2,n);
            end
            

            % Clear score and birth generation
            obj.score = [];
            obj.birthGeneration = [];
        end

        function obj = yShiftKeyPoints(obj, numProfiles)
        % yShiftKeyPoints applies a random vertical shift to each key point
        %   Acts as a function wrapper to obj.single_yShiftKeyPoints
        %   Function wrapper is required to handle calls when the parent
        %   object is an array (multiplicity > 1) or when multiple profiles
        %   need to be generated.
        %
        % obj = obj.yShiftKeyPoints() creates a "y shifted" profile for 
        %   each element in obj.
        %
        % obj = obj.yShiftKeyPoints(numProfiles) creates the specified
        %   number of "y shifted" profiles for each element in obj.
        %
        %   See also single_yShiftKeyPoints
            if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n) = obj(n).single_yShiftKeyPoints();
            end

            if numProfiles > 1
                obj = obj + obj.yShiftKeyPoints(numProfiles-1);
            end
        end

        function obj = single_yShiftKeyPoints(obj)
        % single_yShiftKeyPoints implements y shifting key points 
        %   The "single" version of this function is the actual y shift
        %   implementation for when object multiplicity is 1. Higher 
        %   multiplicity calls are handled by the above wrapper function.
        %
        % See also yShiftKeyPoints

            r = obj.yShiftConstant;

            % Recalculate first key point
            minBound = max(obj.minStartPercentage*obj.maxPressure,obj.keyPoints(2,1) * (1-r));
            maxBound = min(obj.maxStartPercentage*obj.maxPressure,obj.keyPoints(2,1) * (1+r));
            newKeyPoints(1) = obj.randInRange(minBound,maxBound);

            % Recalculate remaining key points
            for n = 2:size(obj.keyPoints,2)-1
                minBound = max(0,obj.keyPoints(2,n) * (1-r));
                maxBound = max(min(obj.maxNextPeak(n),obj.keyPoints(2,n) * (1+r)),obj.maxPressure/6);
                newKeyPoints(n) = obj.randInRange(minBound,maxBound);
            end
            newKeyPoints(end+1) = 0;

            % Resize/stretch coords to new key point locations
            for n = 1:size(obj.keyPoints,2)-1
                % Collect old coordinate data
                oldMin = min(obj.keyPoints(2,n:n+1));
                oldMax = max(obj.keyPoints(2,n:n+1));
                oldRange = oldMax-oldMin;
                oldSlope = obj.keyPoints(2,n+1)-obj.keyPoints(2,n);

                % Collect new coordinate data
                newMin = min(newKeyPoints(n:n+1));
                newMax = max(newKeyPoints(n:n+1));
                newRange = newMax-newMin;
                newSlope = newKeyPoints(n+1)-newKeyPoints(n);

                % Get relevant coord address range
                xRange = obj.keyPoints(1,n):obj.keyPoints(1,n+1)-1;

                % Shift down so lowest coord is zero, flip if needed
                if oldSlope * newSlope > 0
                    obj.coords(xRange) = obj.coords(xRange) - oldMin;
                else
                    obj.coords(xRange) = oldMax - obj.coords(xRange);
                end
                
                % Resize coords
                obj.coords(xRange) = obj.coords(xRange) .* (newRange / oldRange);
                
                % Return coords to new location
                obj.coords(xRange) = obj.coords(xRange) + newMin;
            end

            % Catch edge case when key points formed horizontal line
            obj.coords(find(isnan(obj.coords)|isinf(obj.coords))) = 0;            
            obj.keyPoints(2,:) = newKeyPoints;
            obj = obj.interpolateCoords();
        end

        function obj = xyShiftKeyPoints(obj, numProfiles)
% % %         % yShiftKeyPoints applies a random vertical shift to each key point
% % %         %   Acts as a function wrapper to obj.single_yShiftKeyPoints
% % %         %   Function wrapper is required to handle calls when the parent
% % %         %   object is an array (multiplicity > 1) or when multiple profiles
% % %         %   need to be generated.
% % %         %
% % %         % obj = obj.yShiftKeyPoints() creates a "y shifted" profile for 
% % %         %   each element in obj.
% % %         %
% % %         % obj = obj.yShiftKeyPoints(numProfiles) creates the specified
% % %         %   number of "y shifted" profiles for each element in obj.
% % %         %
% % %         %   See also single_yShiftKeyPoints
            if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n) = obj(n).single_xShiftKeyPoints().single_yShiftKeyPoints();
            end

            if numProfiles > 1
                obj = obj + obj.yShiftKeyPoints(numProfiles-1);
            end
        end

        function obj = changeKeyPointsNum(obj, numProfiles)
        % changeKeyPointsNum changes the number of key points in the curve
        %   Acts as a function wrapper to obj.single_changeKeyPointsNum
        %   Function wrapper is required to handle calls when the parent
        %   object is an array (multiplicity > 1) or when multiple profiles
        %   need to be generated.
        %   Changes number of key points to any number specified in the
        %   array numKeyPointsDomain, excluding the current value.
        %
        % obj = obj.changeKeyPointsNum() creates a curve for each element
        %   in obj with a different number of key points.
        %
        % obj = obj.xShiftKeyPoints(numProfiles) creates the specified
        %   number of profiles for each element in obj, each with a
        %   randomly determined number of key points.
        %
        % See also single_changeKeyPointsNum
            if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n) = obj(n).single_changeKeyPointsNum();
            end

            if numProfiles > 1
                obj = obj + obj.changeKeyPointsNum(numProfiles-1);
            end
        end

        function obj = single_changeKeyPointsNum(obj)
        % single_changeKeyPointsNum implements changing # of key points
        %   The "single" version of this function is the implementation 
        %   for when object multiplicity is 1. Higher multiplicity calls 
        %   are handled by the above wrapper function.
        %
        % See also changeKeyPointsNum

            if numel(obj.numKeyPointsDomain) == 1
                error("numKeyPointsDomain only has 1 element.")
            else
                nKPD = obj.numKeyPointsDomain+2;
                nKPD = nKPD(nKPD ~= size(obj.keyPoints,2));
                numKeyPoints = nKPD(randi([1,numel(nKPD)]));
                obj.keyPoints = round(imresize(obj.keyPoints,[2,numKeyPoints],"bilinear"));
                obj.keyPoints(1,1) = 1;
                obj.keyPoints(1,end) = numel(obj.coords);
                obj.keyPoints(2,end) = 0;
                obj = obj.xShiftKeyPoints();
            end
        end


        function obj = nudge(obj, numProfiles)
            if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n) = obj(n).single_nudge();
            end

            if numProfiles > 1
                obj = obj + obj.nudge(numProfiles-1);
            end
        end

        function obj = single_nudge(obj)
            
            nudgePoint = randi([1,size(obj.keyPoints,2)-1]);
            newKeyPoints = obj.keyPoints(2,:);
            newKeyPoints(nudgePoint) = newKeyPoints(nudgePoint) + (randi([0,1])*2-1)*obj.maxPressure/30;

            % Resize/stretch coords to new key point locations
            for n = max(1,nudgePoint-1):nudgePoint
                % Collect old coordinate data
                oldMin = min(obj.keyPoints(2,n:n+1));
                oldMax = max(obj.keyPoints(2,n:n+1));
                oldRange = oldMax-oldMin;
                oldSlope = obj.keyPoints(2,n+1)-obj.keyPoints(2,n);

                % Collect new coordinate data
                newMin = min(newKeyPoints(n:n+1));
                newMax = max(newKeyPoints(n:n+1));
                newRange = newMax-newMin;
                newSlope = newKeyPoints(n+1)-newKeyPoints(n);

                % Get relevant coord address range
                xRange = obj.keyPoints(1,n):obj.keyPoints(1,n+1)-1;

                % Shift down so lowest coord is zero, flip if needed
                if oldSlope * newSlope > 0
                    obj.coords(xRange) = obj.coords(xRange) - oldMin;
                else
                    obj.coords(xRange) = oldMax - obj.coords(xRange);
                end
                
                % Resize coords
                obj.coords(xRange) = obj.coords(xRange) .* (newRange / oldRange);
                
                % Return coords to new location
                obj.coords(xRange) = obj.coords(xRange) + newMin;
            end

            % Catch edge case when key points formed horizontal line
            obj.coords(find(isnan(obj.coords)|isinf(obj.coords))) = 0;            
            obj.keyPoints(2,:) = newKeyPoints;
            obj = obj.interpolateCoords();
        end


        function obj = reinterpolateCoords(obj,numProfiles)
        % reinterpolateCoords changes the interpolation between key points
        %   Clears the current interpolation between key points and
        %   randomly generates a new one. Works for any multiplicity of obj
        %   and can generate a given number of profiles. It accomplishes
        %   this by resetting the coords vector and then calling
        %   obj.interpolateCoords for each element in obj.
        %
        % obj = obj.reinterpolateCoords() creates a reinterpolated version
        %   of each element in obj.
        %
        % obj = obj.reinterpolateCoords(numProfiles) creates the specified
        %   number of reinterpolated profiles for each element in obj.
        %
        % See also interpolateCoords

            if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n).coords = zeros(1,obj(n).keyPoints(1,end));
                obj(n) = obj(n).interpolateCoords();
            end

            if numProfiles > 1
                obj = obj + obj.reinterpolateCoords(numProfiles-1);
            end
        end
        
        function obj = changeKeyPoints(obj,numProfiles)
           if nargin < 2
                numProfiles = 1;
            end
            for n = 1:numel(obj)
                obj(n) = obj(n).single_changeKeyPoints();
            end

            if numProfiles > 1
                obj = obj + obj.changeKeyPoints(numProfiles-1);
            end
        end

        function obj = single_changeKeyPoints(obj)
            r = obj.xShiftConstant;

            % Anonymous function calculates bounds each key point can shift
            calcCoordBound = @(n,r) round(obj.keyPoints(1,n)*r+obj.keyPoints(1,n+1)*(1-r));
            
            % Calculate shift for all key points (excluding first and last)
            for n = 2:size(obj.keyPoints,2)-1
                % Calculate x shift amount
                newX = randi([calcCoordBound(n-1,r),calcCoordBound(n,1-r)]);
                dx = newX - obj.keyPoints(1,n);
                
                obj.keyPoints(1,n) = newX;
                obj.keyPoints(2,n) = obj.coords(newX);
            end


        end

        function obj = interpolateCoords(obj)
        % interpolateCoords interpolates the coords between each key point
        %   Only works when multiplicity is 1. Only sets values that have 
        %   not been set. To come up with a new interpolation, instead call
        %   obj.reinterpolateCoords
        %
        % See also reinterpolateCoords

            % Set the coords for each key point
            for n = 1:size(obj.keyPoints,2)-1
                obj.coords(obj.keyPoints(1,n)) = obj.keyPoints(2,n);
            end

            % Temporarily set the last coordinate to 1 to mark it as set
            obj.coords(end) = 1;

            % Loop will iterate until the number of set coords (not equal
            % to zero) equals the overall number of coords. In other words,
            % keep setting coords until all coords are set.
            while true
                setCoords = find(obj.coords~=0);
                if numel(setCoords) == numel(obj.coords)
                    break;
                end

                % Anonymous function identifies nonadjacent set coordinates
                % which indicates there are unset coordinates in the middle
                hasGap = @(n) setCoords(n) ~= setCoords(n+1)-1;

                % For loop identifies the first nonadjacent set
                % coordinates, finds the midpoint between them, and
                % randomly sets the y coordinate of that midpoint to
                % somewhere in between the already set coords. This
                % stops the resulting curve from doubling back on itself.
                for n = 1:numel(setCoords)
                    if hasGap(n)
                        gap = round((setCoords(n)+setCoords(n+1))/2);
                        minBound = min(obj.coords(setCoords(n)),obj.coords(setCoords(n+1)));
                        maxBound = max(obj.coords(setCoords(n)),obj.coords(setCoords(n+1)));
                        obj.coords(gap) = obj.randInRange(minBound,maxBound);
                        break;
                    end
                end
            end

            % Reset the ending coord back to zero
            obj.coords(end) = 0;

            % Clear score and birth generation
            obj.score = [];
            obj.birthGeneration = [];
        end

        function mutant = mutate(obj,id,numProfiles)
            if nargin < 3
                numProfiles = 1;
            end
            switch id
                case "nf"
                    mutant = obj.newFrom(numProfiles);
                case "xs"
                    mutant = obj.xShiftKeyPoints(numProfiles);
                case "ys"
                    mutant = obj.yShiftKeyPoints(numProfiles);
                case "xy"
                    mutant = obj.xyShiftKeyPoints(numProfiles);
                case "cn"
                    mutant = obj.changeKeyPointsNum(numProfiles);
                case "nd"
                    mutant = obj.nudge(numProfiles);
                case "rc"
                    mutant = obj.reinterpolateCoords(numProfiles);
            end

        end

        function out = randInRange(obj,minBound,maxBound)
            if obj.yIsInt
                out = randi(round([minBound,max(minBound,maxBound)]));
            else
                out = (maxBound-minBound)*rand() + minBound;
            end
        end

        function iseq = eq(a,b)
            iseq = isequal(a.coords,b.coords);
        end


        function obj = setShowKeyPoints(obj,value)
            for n = 1:numel(obj)
                obj(n).showKeyPoints = value;
            end
        end
        
        function showPlot(obj,lineWidth,opacity)
            if nargin < 2
                lineWidth = ones(size(obj)); 
            end
            if numel(lineWidth) < numel(obj)
                lineWidth = repmat(lineWidth(1),size(obj));
            end
            if nargin < 3
                opacity = ones(size(obj));
            end
            if numel(opacity) < numel(obj)
                opacity = repmat(opacity(1),size(obj));
            end
            
            markerColors = zeros([numel(obj),3]);

            for n = 1:numel(obj)
                plt = plot(obj(n).coords,LineWidth=lineWidth(n));
                plt.Color(4) = opacity(n);
                markerColors(n,:) = plt.Color;
                hold on

            end
            
            for n = 1:numel(obj)
                if obj(n).showKeyPoints
                    s = scatter(obj(n).keyPoints(1,:),obj(n).keyPoints(2,:), ...
                        'filled','MarkerFaceAlpha',opacity(n),['' ...
                        'MarkerFaceColor'],markerColors(n,:));
                    s.SizeData = lineWidth(n) *36;
                end
            end

            xlim([1,numel(obj(1).coords)])
            ylim([0,obj(1).maxPressure])
        end

    end


end
        
