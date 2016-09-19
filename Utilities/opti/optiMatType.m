function matType = optiMatType(mat)
%Determine matrix format (S sym, U triu, L tril, D diag, N no idea)

if(optiIsSymmetric(mat))
    matType = 'S';
    %Check if also diagonal,even if all zero
    if(~isempty(which('isdiag')))
        if(isdiag(mat)) %matlab method
            matType = 'D';
        end
    else
        if(optiIsTriULD(mat) == 'D') %lazy for now
            matType = 'D';
        end
    end
else
    matType = optiIsTriULD(mat);
end

