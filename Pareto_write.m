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


payload = [-1, -2, -3]; %Codes for payload mass
D_1 = [-4, -5, -6]; %Codes for first-stage diameter 
D_2 = [-7, -8]; %Codes for second-stage diameter
alt = [-9, 10, 11]; %Codes for desired altitude
%dv_1 = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]; %Codes for dv % split for stage 1
%coast_lim = [23, 24, 25, 26, 27]; %Codes for coasting upper limit

x = (combvec(payload, D_1, D_2, alt));% dv_1, coast_lim)); % get all possible pareto combinations


writematrix(x, 'combos.txt') %write to file

%% For exporting with zeroes (possibly useful in the future)%% 
% for r = 1:length(x(1,:))
%    for y = 1:length(x(:,1))
%       combos{y,r} = num2str(x(y,r));
%    end
% end
% combos = replace(combos, '-', '0');
% writecell(combos)
%% ---- %%


