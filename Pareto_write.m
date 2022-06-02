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
% ordered list of keycode combinations used in Pareto analysis
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

payload = [-1, -2, -3];
D_1 = [-4, -5, -6];
D_2 = [-7, -8];
alt = [-9, 10, 11];
dv_1 = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
coast_lim = [23, 24, 25, 26, 27];

x = (combvec(payload, D_1, D_2, alt, dv_1, coast_lim));

%% For exporting with zeroes %%
% for r = 1:length(x(1,:))
%    for y = 1:length(x(:,1))
%       combos{y,r} = num2str(x(y,r));
%    end
% end
% combos = replace(combos, '-', '0');
% writecell(combos)
%% ---- %%


writematrix(x, 'combos.txt')



