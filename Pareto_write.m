function Pareto_write
%% Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:
% Pareto_write
%
% Date of Creation:
% 06/22/2022
%
% Author(s):
% Justin Kruse
% Thomas Beimrohr
% 
%
% Description:
% For pareto analysis, codes are used to represent different combinations 
% of design choices and metrics. This function creates a .txt file with an
% ordered list of keycode combinations used in Pareto analysis. "-" is used
% in place of a zero.
%
% Inputs:
% None
%
% Outputs:
% Writes to 'combos.txt'
%
%
% Notes:
% YOU WILL NEED 'DEEP LEARNING TOOLBOX' TO USE THIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

payload = [-1]; %Codes for payload mass
D_1 = [-2,-3,-4]; %Codes for first-stage diameter 
D_2 = [-5,-6,-7]; %Codes for second-stage diameter
alt = [-8,-9,10]; %Codes for desired altitude
dv_1 = [11,12]; %Codes for dv % split for stage 1
coast_lim = [13]; %Codes for coasting upper limit

x = (combvec(payload, D_1, D_2, alt,dv_1,coast_lim))'; % get all possible pareto combinations

for i = 1:length(x)
    combo(i,:) = num2str(x(i,:));
    temp = replace(combo(i,:), '-', '0');
    Combinations(i,:) = temp(find(~isspace(temp)));
end

writematrix(Combinations, 'combos.txt') %write to file
