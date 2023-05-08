function [prob,opts] = QCQP2NLP(prob,opts)
%QCQP2NLP  Converts a LP / BILP / MILP / QCQP / MIQP / MIQCQP to NLP / MINLP
% [prob,opts] = QCQP2NLP(prob,opts)

%   Copyright (C) 2012 Jonathan Currie (Control Engineering)

if(~isfield(prob,'type') || isempty(prob.type))
    error('This function is not for user use, and should only be called from OPTI');
end

ndec = length(prob.f);
nqc = prob.sizes.nqc;

%All require converting objective from LP / QP to NLP
H = prob.H; f = prob.f; bias = prob.objbias; prob.H = []; prob.f = []; prob.objbias = [];
if(isempty(bias)), bias = 0; end
[r,c] = size(H);
if(r+c == 0) %Check for LP
    prob.fun = @(x) f'*x + bias;
    %Gradient
    prob.f = @(x) f;
    %Hessian
    Hc = spalloc(ndec,ndec,0);
    prob.H = @(x,sigma,lambda) Hc;
    prob.Hstr = @() Hc;
else %Must be QP
    prob.fun = @(x) 0.5*x'*H*x + f'*x + bias;
    %Const H Part
    Hc = 0.5*(H + H');
    %Gradient
    prob.f = @(x) (Hc*x + f)'; %transpose to keep derivative checker happy
    %Hessian
    prob.H = @(x,sigma,lambda) sparse(sigma*Hc);
    %Hessian Structure is based on the sum of H and Q
    if(nqc == 1)
        str = H + prob.Q;
    elseif(nqc > 1)
        str = H;
        for i = 1:length(prob.Q)
            str = str + prob.Q{i};
        end
    else
        str = H;
    end
    prob.Hstr = @() sparse(double(str ~= 0));
end

%Only QC requires conversion to nlcon
if(nqc > 0)
    Q = prob.Q; l = prob.l; rl = prob.qrl; ru = prob.qru;
    prob.Q = []; prob.l = []; prob.qrl = []; prob.qru = [];
    %Easy one
    if(nqc == 1) 
        %Const Q Part
        Qc = Q + Q';
        prob.nlcon = @(x) x'*Q*x + l'*x;
        prob.nljac = @(x) (Qc*x + l)';
        prob.nljacstr = @() sparse(ones(1,ndec)); %assume full
        %Update Hessian
        prob.H = @(x,sigma,lambda) sparse(sigma*Hc + lambda(1)*Qc);
    %Harder one
    else
        prob.nlcon = @(x) mQC(x,Q,l,nqc);
        prob.nljac = @(x) mdQC(x,Q,l,nqc,ndec);
        prob.nljacstr = @() sparse(ones(nqc,ndec)); %assume full
        %Update Hessian
        prob.H = @(x,sigma,lambda) mQCHess(sigma,lambda,Hc,Q,nqc);
    end
    %Row Form
    prob.cu = ru;
    prob.cl = rl;
end

%Convert integer arg
if(~isempty(prob.int))
    prob.int = prob.int.str;
end

%Artifically create x0 if it doesn't exist
if(isempty(prob.x0))
    prob.x0 = zeros(size(f));
end

%Rebuild using OPTI constructor
[prob,opts] = opti.buildOpti(prob,opts);



function c = mQC(x,Q,l,nqc)
% Handle to convert cell based QC to nlcon
if(isa(x,'double')), c = zeros(nqc,1); end
for i = 1:nqc
    c(i) = x'*Q{i}*x + l(:,i)'*x;
end

function c = mdQC(x,Q,l,nqc,ndec)
% Handle to convert cell based QC to nljac
c = zeros(nqc,ndec);
for i = 1:nqc
    c(i,:) = ((Q{i} + Q{i}')*x + l(:,i))';
end

function H = mQCHess(sigma,lambda,Hc,Q,nqc)
% Handle to convert cell based QC to enable Hessian Callback
H = sigma*Hc;
for i = 1:nqc
    H = H + lambda(i)*(Q{i} + Q{i}');
end

