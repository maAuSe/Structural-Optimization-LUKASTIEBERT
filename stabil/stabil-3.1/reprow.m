function Matrix=reprow(Matrix,RowSel,nTime,RowInc)

%REPROW   Replicate rows from a matrix.
%
%    Matrix=reprow(Matrix,RowSel,nTime,RowInc)
%   replicates the selected rows from a matrix a number of times and adds them  
%   below the existing rows. The k-th time the increment RowInc is added k times
%   to the copied rows.
%
%   Matrix     Matrix (nRow * nCol)
%   RowSel     Rows to be copied (1 * kRow)
%   nTime      Number of times
%   RowInc     Increments that are added to the different columns of the copied
%              rows (1 * nCol)

% David Dooms
% October 2008

% PREPROCESSING
RowSel=RowSel(:).';
RowInc=RowInc(:).';

nSelRow=length(RowSel);
LastRow=size(Matrix,1);
Matrix=[Matrix; zeros(nTime*nSelRow,size(Matrix,2))];
for k=1:nTime
    Matrix((LastRow+1+(k-1)*nSelRow):(LastRow+k*nSelRow),:)=Matrix(RowSel,:)+k*repmat(RowInc,nSelRow,1);
end
