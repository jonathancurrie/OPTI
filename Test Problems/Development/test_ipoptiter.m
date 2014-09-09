function stop = test_ipoptiter(iter,fval,var)
%Testing IPOPT Iteration callback with extra data

stop = false;
var.x
fprintf('iter: %3d  fval: %12.6g ipr: %g idu: %g\n',iter,fval,var.inf_pr,var.inf_du);