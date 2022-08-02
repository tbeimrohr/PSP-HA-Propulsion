function [newvector] = ModifySize(v,n)

%% Header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: 
% ModifySize
%
% Date of Creation: 
% 03/15/2022
%
% Author(s): 
% Thomas Beimrohr
%
% Description: 
% The purpose of this code is to resize a vector without
% changing the elements of the vector. For example a vector x = [1 2 3]
% can be modified to be x = [1 1.5 2 2.5 3]
%
% Inputs:
% (required) v = original numeric vector
% (required) n = desired size to convert the vector to
%
% Outputs:
% (required) newvector = a vector which contains the desired amount of
% elements and contains the same bounds as the original vector
%
% Notes:
% No known bugs need fixed or improvements need to be made
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of Code

original_length = length(v);
new_rate = ((n - 1)/(original_length - 1))^-1;
newvector = interp1(1:original_length,v,1:new_rate:original_length);