function tri = optiIsTriULD(H)
%Determine if a matrix is tril (L), triu (U), diagonal (D) or neither (N)

%Matlab routines if available
if(~isempty(which('istril')))
    if(istril(H))
        tri = 'L';
    elseif(istriu(H))
        tri = 'U';
    elseif(isdiag(H))
        tri = 'D';
    else
        tri = 'N';
    end
else %JC method
    nnzU = nnz(triu(H,1));
    nnzL = nnz(tril(H,-1));

    if(nnzU == 0 && nnzL > 0)
        tri = 'L';
    elseif(nnzL == 0 && nnzU > 0)
        tri = 'U';
    elseif(nnzU == 0 && nnzL == 0 && nnz(H) > 0)
        tri = 'D';
    else
        tri = 'N';
    end
end
