function rnge=xlsAddr(row,col)
% Build Excel cell address from row, column
%
% RNGE=XLSADDR(ROW,COL) will return an Excel cell address 
% formed from the input ROW,COL values.  Either input may be
% string or numeric and will be converted to canonical form

if isnumeric(col)
  if ~isscalar(col),  error('Input Column Not Scalar'), end
  N=26+(col==26);     % kludge to fix up rollover Z-->AA 
  rnge=num2str('A'+[fix(col/N) rem(col,N)]-1,'%c%c');
  rnge(rnge=='@')=[];   % cleanup for single character
else
  rnge=[col];
end
if isnumeric(row)
  if ~isscalar(row),  error('Input Row Not Scalar'), end
  rnge=[rnge num2str(row,'%d')];
else
  row=num2str(row,'%d');
  if ~all(ismember(row,'0':'9')),  error('Invalid Excel Address: Row not numeric'), end
  rnge=[rnge row];
end