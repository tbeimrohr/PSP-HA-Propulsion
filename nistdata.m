function data = nistdata(species,T,p)
%NISTDATA:  Create tables of thermophysical properties for gases
% Data are read from the NIST chemistry webbook: 
%     http://webbook.nist.gov/chemistry/fluid/
% Input:
%     species: Chemical symbol (e.g. 'H2')
%     T      : Temperature array (K).  Note: must be equally spaced 
%     p      : Pressure array (bar). Must be equally spaced if T is scalar
% Output:
%     data   : Struct with the following fields: 
%       Single values:
%         species : Chemical symbol (e.g. 'H2')
%         Tc      : Critical temperature (K)
%         Pc      : Critical pressure (Pa)
%         Mw      : Molar mass (kg/kmol)
%       Arrays:
%         T       : Temperature (K)
%         P       : Pressure (Pa)
%         Rho     : Density (kmol/m3)
%         V       : Volume (m3/kmol)
%         U       : Internal energy (J/kmol)
%         H       : Enthalpy (J/kmol)
%         S       : Entropy (J/kmol/K)
%         Cv      : Heat capacity at constant volume (J/kmol/k)
%         Cp      : Heat capacity at constant pressure (J/kmol/k)
%         C       : Speed of sound (m/s)
%         JT      : Joule-Thompson coefficient (K/Pa)
%         mu      : Dynamic viscosity (Pa s)
%         k       : Thermal conductivity (W/m/K)
%
% Available species:
% H2, He, N2, O2, Ar, CO, CO2, NH3, H2O
%
% The code is tested mainly for single-phase gas regions, but will 
% provide data for liquid in the liquid region.
%
% In order to avoid unnecessary load on the NIST servers, please save and 
% reuse the data struct, rather than repeatedly reading the same data.
%
% Example:
%  data = nistdata('N2',200:20:400,[1,50:50:300]);
%  surf(data.P*1e-5,data.T,data.C); title('Speed of sound in nitrogen')
%      xlabel('bar');ylabel('K'),zlabel('m/s')
% To find C at 234 bar and 317 K:
%  C = interp2(data.P,data.T,data.C,234e5,317)

    if any(p>1000)
        error('Pressures should not exceed 1000 bar')
    end
    % Species data
    sp = {
    % name, CAS no. ,    Tc ,  Pc     ,acentric,   Mw}        
     'H2O','7732185', 647.  ,220.64e5 ,  -0.344, 18.0153
     'H2' ,'1333740',  33.18, 13.00e5 ,  -0.220,2.01588
     'He' ,'7440597',   5.2 ,  2.274e5,  -0.387,4.002602       
     'Ar' ,'7440371', 150.65,  4.898e5,  -0.004,39.948
     'N2' ,'7727379', 126.19, 33.978e5,  -0.040,  28.013
     'O2' ,'7782447', 154.58, 50.43e5 ,   0.021,31.998
     'CO' ,'630080' , 134.45, 35.0e5  ,   0.049,28.001
     'CO2','124389' , 304.18, 73.8e5  ,   0.099,44.0095       
     'NH3','7664417', 405.4 ,113.0e5  ,   0.250,17.0305
    };
    index = find(strcmp(sp(:,1),species));
    if isempty(index)
        error('Cannot find CAS number for %s',species)
    end
    casno = sp{index,2};
    data.species = species;
    data.Tc = sp{index,3};
    data.pc = sp{index,4};
    data.omega = sp{index,5}; % acentric factor
    data.Mw = sp{index,6};
    
    nT = length(T);
    np = length(p);
    if nT > 1
        dT = T(2)-T(1);
        if ~(all(diff(T)==dT))
            error('The temperature array must be evenly spaced.')
        end        
        type = 'IsoBar';
        isobar = true;
        dT = T(2)-T(1);
        if ~(all(diff(T)==dT))
            error('The temperature array must be evenly spaced.')
        end        
        tstr = sprintf('&TLow=%g&THigh=%g&TInc=%g',T(1),T(end),dT);
    elseif np > 1
        tstr = sprintf('&T=%g',T(1));
        dp = p(2)-p(1);
        if (all(diff(p)==dp))        
            type = 'IsoTherm';
            isobar = false;        
            pstr = sprintf('&PLow=%g&PHigh=%g&PInc=%g',p(1),p(end),dp);
        else
            error('If T is scalar P must be equally spaced')
        end
    else
        error('T and P cannot both be scalars')        
    end
    str1 = 'http://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C';
    str2 = [casno,'&Type=',type,'&Digits=5'];
    str3 = '&RefState=DEF&TUnit=K&PUnit=bar&DUnit=mol%2Fl&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm';    
    fields = {'T','P','Rho','V','U','H','S','Cv','Cp','C','JT','mu','k'}; 
    % Scale factors to SI values (using kmol instead of mol):
    scale  = [ 1 ,1e5,  1  , 1 ,1e6,1e6,1e3, 1e3, 1e3, 1 ,1e-5,1,1];
    nf = length(fields);
    if isobar
        for j = 1:np
            pstr = sprintf('&P=%g',p(j));
            url = [str1,str2,pstr,tstr,str3];
            table = parse_nist(url);
            % NIST inserts two extra rows showing liquid and vapour
            % properties at the boiling temperature.  We remove these 
            % to keep the number of rows fixed:
            deltaT = diff(table(:,1));
            twophase = find(deltaT==0,1);
            if ~isempty(twophase)
                table(twophase+[0,1],:) = [];
            end
            for k = 1:nf
                data.(fields{k})(:,j) = table(:,k)*scale(k);
            end
        end
    else  % Isotherm
        url = [str1,str2,pstr,tstr,str3];
        table = parse_nist(url);
        % When the isotherm range spans the critical pressure, NIST adds
        % two rows for vapour and supercritical. We remove these 
        % superfluous rows:
        deltaP = diff(table(:,2));
        repeats = find(deltaP==0);
        for j = length(repeats):-1:1
            table(repeats(j),:) = [];
        end  
        for k = 1:nf
            data.(fields{k}) = table(:,k)*scale(k);
        end         
    end    
end


function txt = read_nist(url)
% Read and interpret data from the NIST fluids webbook
    s = regexp(version,'\.','split');
    version_number = str2double([s{1},s{2}]);
    if version_number >=93  
        % Version R2017b or newer should be able to read secure urls 
        % (https), but my Matlab version is too old to test this.  
		    % Thanks to Martijn van Heumen for his comment with this solution:
        options = weboptions('RequestMethod','get','ArrayFormat', ...
          'csv','ContentType','text'); 
        txt = webread(url,options);
    else % Use Python to read the url   
        if ~exist('read_url.py','file')
            fprintf(no_read_url)
            error('read_url.py not found')
        end        
        try
            if count(py.sys.path,'') == 0
                insert(py.sys.path,int32(0),'');
            end
            txt = char(py.read_url.read_url(url));
            % Python 3 returns linefeed and tab as literal '\n' and '\t'
            % Replcae by linefeed and tab characters:
            txt = regexprep(txt,'\\n','\n');
            txt = regexprep(txt,'\\t','\t');
        catch
            fprintf(no_python);
            error('python not found')
        end
    end
end

function table = parse_nist(url)
    txt = read_nist(url);
    txt = regexp(txt,'\n','split');
    nlines = length(txt)-1;
    tabpos = find(txt{2}==9);
    ncols = length(tabpos)+1;

    % Read data into array table
    table = zeros(nlines-1,ncols)*NaN;
    for i = 2:nlines
        if isempty(txt{i})
            continue
        end
        tabpos = find(txt{i}==9);
        colstart = [1,tabpos+1];
        colend = [tabpos-1,length(txt{i})];
        for j = 1:ncols
            table(i-1,j) = str2double(txt{i}(colstart(j):colend(j)));
        end
    end 
end    

function t = no_read_url
t = ['\nERROR: Could not find read_url.py\n\n',...
    'Matlab versions older than R2017b cannot read secure urls ',...
    '(httpps),\nso I use Python for this task. Unfortunately, '...
    'this seems to \nhave failed. Note that\n',...
    '   1) Python (2.7 or 3) must be installed.\n   2) The file ',...
    'read_url.py must exist in the Matlab path.\n      This file is ',...
    ' distributed along with nistdata.m\n',...
    'You may have to restart Matlab after fixing this problem.\n\n'];
end

function t = no_python
t = ['\nERROR: Could not find python\n\n',...
    'My Matlab version (2014b) cannot read secure urls (httpps), so ',...
    'I use Python \nfor this task. Unfortunately, this seems to have ',...
    'failed in this case.\nNote that Python (2.7 or 3) must be ',...
    'installed.\n\nIf you use Windows, the environment variable ',...
    '"path" must contain the folder \nfor python.exe. ',...
    'You may have to restart Matlab after fixing this problem.\n\n'];
end
