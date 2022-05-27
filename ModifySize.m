function [newvector] = ModifySize(v,n)
original_length = length(v);
new_rate = ((n - 1)/(original_length - 1))^-1;
newvector = interp1(1:original_length,v,1:new_rate:original_length);