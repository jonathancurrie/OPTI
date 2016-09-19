function issym = optiIsSymmetric(H)
%Returns true if a matrix is symmetric, using MATLAB routine if available

%Check if MATLAB version is available (>= R2014a)
if(~isempty(which('issymetric')))
    issym = issymmetric(H);
else %jc method
    if(nnz(triu(H,1) - tril(H,-1).') == 0)
        issym = true;
    else
        issym = false;
    end
end