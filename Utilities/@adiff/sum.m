function c = sum(a,dim)
% SUM of the adiff object. 

c = adiff( sum(a.x), sum(a.dx), a.root);