function displayBuilder(B)
%Display Builder Problem Information
%
%   Called By SYMBUILDER Class

%   Copyright (C) 2012 Jonathan Currie (Control Engineering)

disp('------------------------------------------------------');
disp('SymBuilder Object');

%Display Based on Build Status
switch(B.bldstat)
    case 'unbuilt'
        fprintf(' UNBUILT with:\n');
        fprintf(' - %2d objective(s)\n',B.noObjs);
        fprintf(' - %2d constraint(s)\n',B.noCons);
        fprintf(' - %2d bound declaration(s)\n',size(B.bnds,1));
        fprintf(' - %2d integer declaration(s)\n',size(B.vartypes,1));
        fprintf(' - %2d expression(s)\n',size(B.exprsn,1));
        fprintf(' - %2d constant(s)\n',size(B.constnt,1));       
    case {'draft'}        
        fprintf(' DRAFT in %1.3fs with:\n',B.bldtime);
    case 'built'
        fprintf(' BUILT in %1.3fs with:\n',B.bldtime);
end

%If object is built, display properties
if(~strcmp(B.bldstat,'unbuilt'))
    fprintf(' - %2d variables\n',length(B.vars));
    if(B.noObjs)
        fprintf(' - %2d objective\n',B.noObjs);
        fprintf('      - %2d linear\n',sum(B.objlin == 1));
        if(~isempty(B.hess))
            fprintf('      - %2d quadratic\n',sum(B.objlin == 2));
        end
        fprintf('      - %2d nonlinear\n',sum(B.objlin == 3));
    else
        fprintf(' - No objective\n');
    end
    if(B.noCons)
        fprintf(' - %2d constraint(s)\n',B.noCons);
        fprintf('      - %2d linear\n',sum(B.conlin == 1));
        if(~isempty(B.hess))
            fprintf('      - %2d quadratic\n',sum(B.conlin == 2));
        end
        fprintf('      - %2d nonlinear\n',sum(B.conlin == 3));
    else
        fprintf(' - No constraints\n');
    end
    fprintf(' - %2d bound(s)\n',sum(~isinf(B.lb)) + sum(~isinf(B.ub)));
    int = sum(B.xtype ~= 'C');
    fprintf(' - %2d integer variable(s)\n',int);
    if(int)
        fprintf('      - %2d integer\n',sum(B.xtype == 'I'));
        fprintf('      - %2d binary\n',sum(B.xtype == 'B'));
    end
    
%     if(any(B.objlin == 1))
%         fprintf('\n - Linear Objective:\n');
%         lin = B.sobj(B.objlin == 1);
%         for i = 1:length(lin)
%            fprintf('%40s\n',char(lin(i)));
%         end
%     end
%     if(any(B.objlin == 2))
%         fprintf('\n - Quadratic Objective:\n');
%         lin = B.sobj(B.objlin == 2);
%         for i = 1:length(lin)
%            fprintf('%40s\n',char(lin(i)));
%         end
%     end
%     if(any(B.objlin == 3))
%         fprintf('\n - Nonlinear Objective:\n');
%         nl = B.sobj(B.objlin == 3);
%         for i = 1:length(nl)
%             fprintf('%40s\n',char(nl(i)));
%         end
%     end
% 
%     if(any(B.conlin == 1))
%         fprintf('\n - Linear Constraints [%d nz]:\n',nnz(sparse(double(B.jac(B.conlin==1,:)))));
%         ind = B.conlin == 1;
%         lin = B.sobj(ind);
%         lwr = B.cl(ind);
%         upr = B.cu(ind);
%         for i = 1:length(lwr)
%             [s,r] = getConType(lwr(i),upr(i));                
%             fprintf('%40s %2s %1.5g\n',char(lin(i)),s,r);
%         end
%     end
%     if(any(B.conlin == 2))
%         fprintf('\n - Quadratic Constraints:\n');
%         ind = B.conlin == 2;
%         q = B.sobj(ind);
%         lwr = B.cl(ind);
%         upr = B.cu(ind);
%         for i = 1:length(lwr)
%             [s,r] = getConType(lwr(i),upr(i));                
%             fprintf('%40s %2s %1.5g\n',char(q(i)),s,r);
%         end
%     end
%     if(any(B.conlin == 3))
%         fprintf('\n - Nonlinear Constraints:\n');
%         ind = B.conlin == 3;
%         nl = B.sobj(ind);
%         lwr = B.cl(ind);
%         upr = B.cu(ind);
%         for i = 1:length(lwr)
%             [s,r] = getConType(lwr(i),upr(i));                
%             fprintf('%40s %2s %1.5g\n',char(nl(i)),s,r);
%         end
%     end
end
disp('------------------------------------------------------');


function [s,r] = getConType(cl,cu)
%Determine constraint type based on supplied args
if(isinf(cl))
    s = '<='; r = cu;
elseif(isinf(cu))
    s = '>='; r = cl;
elseif(cl == cu)
    s = '=='; r = cl;
else
    error('This interface does not yet support double constraints');
end
        
