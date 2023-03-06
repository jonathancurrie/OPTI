function p = optiWriteSDPA(filename,f,A,b,lb,ub,sdcone,xtype,dense)
%OPTIWRITESDPA  Write an SDPA (.dat-s or .dat) file from Matlab Values
%
% optiWriteSDPA(filename,f,A,b,lb,ub,sdcone) writes the sparse SDPA file 
% specified by filename and f (objective), A,b, (linear inequalities), lb,ub
% (bounds) and  sdcone (semidefinite constraints). If you do not specify a 
% file a extension, it will default to .dat-s.
%
% optiWriteSDPA(filename,...,xtype) allows the user to specify 
% variable integrality ('C','I','B'). The default is all continuous.
%
% optiWriteSDPA(filename,...,dense) allows the user to specify 
% whether the file is dense or not. dense = 1 indicates a dense SDPA file, 
% and the file extensions will default to .dat if not provided. The default 
% is sparse.
% 
% *sdcone is a cell array, where each cell is a sparse matrix of the
% form [C(:) A0(:) A1(:) ... ], containing symmetric matrices.

%   Copyright (C) 2013 Jonathan Currie (Control Engineering)

if(nargin < 7), error('You must supply at least 7 arguments to this function'); end
if(nargin < 8), dense = 0; end

% Make sure we have file extension!
if(isempty(strfind(filename,'.')))
    filename = [filename '.dat-s'];
end
if(~isempty(sdcone) && ~iscell(sdcone))
    sdcone = {sdcone};
end

if(isempty(f)), error('You must supply an objective vector to this function'); end   
ndec = length(f);
if(~isempty(A))
    if(isempty(b)), error('You must supply both A and b for inequality constraints'); end
    if(size(A,1) ~= length(b)), error('Dimensions of A and b are incorrect'); end
    if(size(A,2) ~= ndec), error('Dimensions of A and f are incorrect'); end
    cA = 1; 
else
    cA = 0; 
end
if(~isempty(lb))
    if(length(lb) ~= ndec), error('f and lb are not the same length'); end
    clb = 1; 
else
    clb = 0; 
end
if(~isempty(ub))
    if(length(ub) ~= ndec), error('f and ub are not the same length'); end
    cub = 1; 
else
    cub = 0; 
end
if(~isempty(b) && any(isinf(b))), error('All elements in b must be finite!'); end

%Check integer constraints
if(~isempty(xtype))
    xtype = lower(xtype);
    if(any(xtype ~= 'c'))
        cint = 1;
        intidx = find(xtype ~= 'c');
        %check binary variables, if so, we must bound
        if(any(xtype == 'b'))
            if(clb==0), lb = -Inf(ndec,1); end
            if(cub==0), ub = Inf(ndec,1); end
            bind = find(xtype == 'b');
            for i = 1:length(bind)
                if(lb(bind(i)) < 0)
                    lb(bind(i)) = 0;
                end
                if(ub(bind(i)) > 1)
                    ub(bind(i)) = 1;
                end
            end
        end 
    else
        cint = 0;
    end
else
    cint = 0;
end

%Open File
fid = fopen(filename,'w+');
if(fid == -1)
    error('Error opening file %s',filename);
end

%Begin writing data
try
    %Write Heading
    if(dense)
        fprintf(fid,'"OPTI SDP Problem [Dense] generated on %s"\n',datestr(now));
    else
        fprintf(fid,'"OPTI SDP Problem [Sparse] generated on %s"\n',datestr(now));
    end
    %Generate block struct
    blockstruct = zeros(1,cA+clb+cub+length(sdcone)); offset=0;
    if(~isempty(lb))
        blockstruct(1) = -sum(~isinf(lb)); offset = offset + 1;
    end
    if(~isempty(ub))
        blockstruct(1+offset) = -sum(~isinf(ub)); offset = offset + 1;
    end
    if(~isempty(A))
        blockstruct(1+offset) = -size(A,1); offset = offset + 1;
    end
    for i = 1:length(sdcone)
        dim = sqrt(size(sdcone{i}(:,1),1));
        if(floor(dim) ~= dim), error('SDCone{%d} does not contain the correct number of rows to form a square matrix',i); end
        if(size(sdcone{i},2) ~= ndec+1), error('SDCone{%d} does not contain the correct number of columns, expected %d',i,ndec+1); end
        blockstruct(i+offset) = dim;
    end
        
    %Write Sizes
    fprintf(fid,'%d\n%d\n',ndec,length(blockstruct));   
    writeVec(fid,blockstruct);
    
    %Write f
    writeVec(fid,-f);
    
    %Write C Matrices
    writeC(fid,lb,ub,b,sdcone,dense);
    %Write Each A Matrix
    for i = 1:ndec
        writeA(fid,i,ndec,lb,ub,A,sdcone,dense);
    end
    %Write Integer Declaration
    if(cint)
        fprintf(fid,'*INTEGER\n');
        for i = 1:length(intidx)
            fprintf(fid,'*%d\n',intidx(i));
        end
    end
    
catch ME
    fclose(fid);
    rethrow(ME);
end
%Close file
fclose(fid);
p = which(filename);


%Write all C matrix blocks [lb,ub,lin,sdcone]
function writeC(fid,lb,ub,b,sdcone,dense)
%Write bounds and linear constraints as diagonal entries
con = 0;
block = 1;
if(dense), fprintf(fid,'{\n'); end
if(~isempty(lb)), writeVec(fid,lb(~isinf(lb)),con,block,dense,1); block = block + 1; end
if(~isempty(ub)), writeVec(fid,-ub(~isinf(ub)),con,block,dense,1); block = block + 1; end
if(~isempty(b)), writeVec(fid,-b,con,block,dense,1); block = block + 1; end
%For each semidefinite block
for i = 1:length(sdcone)
    writeMat(fid,sdcone{i}(:,1),con,i+block-1,dense);        
end
if(dense), fprintf(fid,'}\n'); end
    
%Write all A_i matrix blocks [lb,ub,lin,sdcone]
function writeA(fid,con,ndec,lb,ub,A,sdcone,dense)
block = 1;
if(dense), fprintf(fid,'{\n'); end
if(~isempty(lb) && ~isinf(lb(con)))
    ninf = sum(isinf(lb));
    %Check if this bound is inf
    if(isinf(lb(con)))
        %Generate empty entry
        vec = sparse([],[],[],ndec-ninf,1,0);
    else
        %Generate linear entry for the bound
        nlinf = sum(isinf(lb(1:con)));
        vec = sparse(con-nlinf,1,-1,ndec-ninf,1,1);
    end
    writeVec(fid,vec,con,block,dense,1); 
    block = block + 1; 
end
if(~isempty(ub))
    ninf = sum(isinf(ub));
    %Check if this bound is inf
    if(isinf(ub(con)))
        %Generate empty entry
        vec = sparse([],[],[],ndec-ninf,1,0);
    else
        %Generate linear entry for the bound
        nlinf = sum(isinf(ub(1:con)));
        vec = sparse(con-nlinf,1,1,ndec-ninf,1,1);
    end
    writeVec(fid,vec,con,block,dense,1); 
    block = block + 1; 
end
if(~isempty(A))
    writeVec(fid,A(:,con),con,block,dense,1); 
    block = block + 1; 
end
%For each semidefinite block
for i = 1:length(sdcone)
    writeMat(fid,-sdcone{i}(:,con+1),con,i+block-1,dense);        
end
if(dense), fprintf(fid,'}\n'); end


%Low level matrix/vector writing functions
function writeMat(fid,mat,con,block,dense)
%mat is supplied as a column vector in column major
dim = sqrt(length(mat));
%If dense, write each row
if(dense)
    for i = 1:dim
        ind = i:dim:dim^2;
        writeVec(fid,mat(ind),con,block,dense);
    end
%If sparse, only writing upper triangular elements    
else
    writeSpTriuMat(fid,mat,con,block);
end

function writeVec(fid,vec,con,block,dense,diag)
if(nargin < 3), con = 0; end
if(nargin < 4), block = 0; end
if(nargin < 5), dense = 1; end
if(nargin < 6), diag = 0; end
fmt = '%1.16g';
if(dense)
    %Don't worry about diagonal for dense vectors, as only a vector
    for i = 1:length(vec)-1
        if(~isinf(vec(i)))
            fprintf(fid,[fmt ' '],full(vec(i)));
        end
    end
    if(~isinf(vec(end)))
        fprintf(fid,[fmt '\n'],full(vec(end)));
    else
        fprintf(fid,'\n');
    end
else
    %Ensure a column vector
    if(size(vec,2) > size(vec,1)), vec = vec'; end
    %Find nonzero elements
    [ii,jj,val] = find(vec);
    %Write in each element
    for i = 1:length(ii)
        if(diag)
            fprintf(fid,['%d %d %d %d ',fmt,'\n'],con,block,ii(i),ii(i),val(i));
        else
            fprintf(fid,['%d %d %d %d ',fmt,'\n'],con,block,ii(i),jj(i),val(i));
        end
    end    
end

function writeSpTriuMat(fid,mat,con,block)
%mat is supplied as a sparse column vector, all elements, column major
fmt = '%1.16g';
dim = sqrt(length(mat));
[ii,~,val] = find(mat);
%Determine row, column from row index in column
r = mod(ii-1,dim)+1;
c = (ii - r)/dim +1;
for i = 1:length(ii)    
    %Write in only upper tri entries
    if(r(i) <= c(i))
        fprintf(fid,['%d %d %d %d ',fmt,'\n'],con,block,r(i),c(i),val(i));
    end
end


