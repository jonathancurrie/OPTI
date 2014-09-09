function x0 = opti_guessX0(prob,ndec)
%Try choose a better x0 than zeros() if we have all finite bounds

if(~isempty(prob.lb) && ~isempty(prob.ub) && all(~isinf(prob.lb)) && all(~isinf(prob.ub)))
    x0 = (prob.ub-prob.lb)./2;
else
    x0 = zeros(ndec,1);
end