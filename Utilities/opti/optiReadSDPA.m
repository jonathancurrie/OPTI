function [f,A,b,sdcone,xtype] = optiReadSDPA(filename,dense,print)
%OPTIREADSDPA  Read an SDPA (.dat-s or .dat) file and convert to Matlab Values
%
% [f,A,b,sdcone] = optiReadSDPA(filename) reads the sparse file specified by 
% filename and return the semidefinite problem in OPTI format. If you do 
% not specify a file a extension, it will default to .dat-s.
%
% [f,A,b,sdcone] = optReadSDPA(filename,dense) allows the user to specify 
% whether the file is dense or not. dense = 1 indicates a dense SDPA file, 
% and the file extensions will default to .dat if not provided. The default 
% is sparse.
%
% [f,A,b,sdcone] = optiReadSDPA(filename,dense,print) allows print out of
% reading progress.
%
% [f,A,b,sdcone,xtype] = optiReadSDPA(..) returns variable integrality.
% 
% *sdcone will be a cell array, where each cell is a sparse matrix of the
% form [C(:) A0(:) A1(:) ... ], containing symmetric matrices.

%   Copyright (C) 2013 Jonathan Currie (I2C2)
%#ok<*AGROW>

if(nargin < 2), dense = 0; end
if(nargin < 3), print = 0; end

% Make sure we have file extension!
if(isempty(strfind(filename,'.')))
    filename = [filename '.dat-s'];
end

%Open File
fid = fopen(filename);
if(fid == -1)
    error('Error opening file %s',filename);
end

%States 
sNCON = 1;
sNBLOCKS = 2;
sSIZES = 3;
sRHS = 4;
sDNDATA = 5;
sSPDATA = 6;
sINTEGER = 7;

%Default Internals
eof = 0;
state = sNCON; cind = 0; bind = 1; rind = 1;
f = []; A = []; b = []; sdcone = []; xtype = [];

if(print)
    ind = strfind(filename,filesep);
    if(isempty(ind))
        file = filename;
    else
        file = filename(ind(end)+1:end);
    end
    if(dense)
        fprintf('Reading Dense SDPA File ''%s''...\n',file);
    else
        fprintf('Reading Sparse SDPA File ''%s''...\n',file);
    end
end
try
    %Main Loop
    while(~eof)
        %Read DAT File String
        str = fgetl(fid); 
        %Check for empty
        if(~ischar(str) && str == -1)
            eof = 1;
            continue;
        end
        %Check for comment, skip if found
        if(state ~= sINTEGER && (str(1) == '*' || str(1) == '"'))
            if(~isempty(strfind(str,'*INTEGER')))
                state = sINTEGER;
            end
            continue;
        end
        %Check for equals (Sonja files)
        idx = strfind(str,'=');
        if(~isempty(idx))
            str = str(1:idx-1);
        end
        %Replace { } ( ) [ ] , ; : with spaces
        try
            str = regexprep(str,{'[',']',',','{','}','(',')'},' ');
        catch
            error('Error in SDPA File with "%s"',str);
        end
        %Decode based on state
        switch(state)
            case sNCON
                M = sscanf(str,'%d',1);
                xtype = repmat('C',1,M);
                %Check for errors
                if(M < 1 || M > 1e9), error('The number of constraints is negative or > 1e9'); end
                %Print what we found
                if(print), fprintf('Variables:              %6d\n',M); end 
                state = sNBLOCKS;
            case sNBLOCKS
                N = sscanf(str,'%d',1);
                %Check for errors
                if(N < 1 || N > 1e9), error('The number of blocks is negative or > 1e9'); end
                %Print what we found
                if(print), fprintf('Number of Blocks:       %6d\n',N); end 
                state = sSIZES;
            case sSIZES
                blocksize = sscanf(str,'%d');
                %Check length = number of blocks
                if(length(blocksize) ~= N), error('The number of blocks does not match the number of entries of block sizes'); end                
                if(any(blocksize == 0) || any(blocksize > 1e9)), error('A block size is 0 or > 1e9'); end
                %Print Total Dims
                if(print)
                    fprintf('  Linear:               %6d\n',length(blocksize(blocksize < 0)));
                    fprintf('  Semidefinite:         %6d\n',length(blocksize(blocksize > 0)));                    
                    fprintf('Total Block Dimension:  %6d\n',sum(abs(blocksize)));                     
                end
                %Allocate problem data
                sdcone = cell(sum(blocksize > 0),1); j = 1;
                for i = 1:N
                    if(blocksize(i) > 0)
                        if(dense)
                            sdcone{j} = zeros(blocksize(i)^2,M+1); 
                        else
                            sdcone{j} = sparse(blocksize(i)^2,M+1); %not the greatest...
                        end
                        j = j + 1;
                    end
                end                
                %Create linear part
                ind = blocksize < 0;
                if(any(ind))
                    if(dense)
                        A = zeros(sum(abs(blocksize(ind))),M);
                    else
                        A = sparse(sum(abs(blocksize(ind))),M);
                    end
                    b = zeros(sum(abs(blocksize(ind))),1);
                end
                %Create index vector
                blockindex = zeros(N,1);  lp = 1; sdp = 1;
                for i = 1:N
                    if(blocksize(i) < 0)
                        blockindex(i) = lp; lp  = lp + 1;
                    else
                        blockindex(i) = sdp; sdp = sdp + 1;
                    end
                end
                state = sRHS;
            case sRHS
                f = sscanf(str,'%f'); f = -f;
                %Check length = number of constraints
                if(length(f) ~= M), error('The number of RHS elements does not match the number of constraints'); end                
                if(dense)
                    state = sDNDATA;
                else
                    state = sSPDATA;
                end
            %Sparse Data File
            case sSPDATA
                data = sscanf(str,'%d %d %d %d %f');
                if(~isempty(data)) %could be {}
                    %Check for C matrix
                    if(data(1) == 0)
                        %Check block for linear 
                        if(blocksize(data(2)) < 0)
                            if(data(3) ~= data(4)), error('A diagonal entry in C is not on the diagonal!'); end
                            %Find start index by determining sizes of previous LP blocks encountered
                            bsizes = blocksize(1:data(2)-1);
                            stIndex = abs(sum(bsizes(bsizes < 0)));
                            b(stIndex + data(3)) = -data(5);
                        else %semidef
                            sdcone{blockindex(data(2))}(data(3) + (data(4) - 1)*blocksize(data(2)),1) = data(5); 
                            %if off-diagonal, place in transposed position (build symmetric matrices)
                            if(data(3) ~= data(4))
                                sdcone{blockindex(data(2))}(data(4) + (data(3) - 1)*blocksize(data(2)),1) = data(5);
                            end
                        end
                    else
                        %Check block for linear
                        if(blocksize(data(2)) < 0)
                            if(data(3) ~= data(4)), error('A diagonal entry in A is not on the diagonal!'); end
                            %Find start index by determining sizes of previous LP blocks encountered
                            bsizes = blocksize(1:data(2)-1);
                            stIndex = abs(sum(bsizes(bsizes < 0)));
                            A(stIndex + data(3),data(1)) = data(5);
                        else %semidef
                            sdcone{blockindex(data(2))}(data(3) + (data(4) - 1)*blocksize(data(2)),data(1)+1) = -data(5);
                            %if off-diagonal, place in transposed position (build symmetric matrices)
                            if(data(3) ~= data(4))
                                sdcone{blockindex(data(2))}(data(4) + (data(3) - 1)*blocksize(data(2)),data(1)+1) = -data(5);
                            end
                        end
                    end
                end
            %Dense Data File
            case sDNDATA
                data = sscanf(str,'%f');
                if(~isempty(data)) %could be {}
                    %Check for C matrix
                    if(cind == 0)
                        %Check block for linear
                        if(blocksize(bind) < 0)
                            if(length(data) ~= abs(blocksize(bind))), error('Data read in dense entry C (diagonal) is not the correct length!'); end
                            %Find start index by determining sizes of previous LP blocks encountered
                            bsizes = blocksize(1:bind-1);
                            stIndex = abs(sum(bsizes(bsizes < 0)))+1;
                            b(stIndex:stIndex+abs(blocksize(bind))-1) = -data;                        
                        else
                            if(length(data) == blocksize(bind))
                                sdcone{blockindex(bind)}(rind:blocksize(bind):blocksize(bind)^2) = data;
                                rind = rind + 1; 
                            elseif(length(data) == blocksize(bind)^2)
                                %Create matrix and transpose it
                                temp = reshape(data,blocksize(bind),blocksize(bind))';
                                sdcone{blockindex(bind)}(:,1) = temp(:);
                                rind = rind + blocksize(bind);
                            else
                                error('Data read in dense entry C matrix is not the correct length. Expected first row of a matrix, or the entire matrix (all rows and columns, row major).'); 
                            end                          
                        end
                    else
                        %Check block for linear
                        if(blocksize(bind) < 0)
                            if(length(data) ~= abs(blocksize(bind))), error('Data read in dense entry A%d (diagonal) is not the correct length!',cind); end
                            %Find start index by determining sizes of previous LP blocks encountered
                            bsizes = blocksize(1:bind-1);
                            stIndex = abs(sum(bsizes(bsizes < 0)))+1;                            
                            A(stIndex:stIndex+abs(blocksize(bind))-1,cind) = data;
                        else
                            if(length(data) == blocksize(bind))
                                sdcone{blockindex(bind)}(rind + (cind)*blocksize(bind)^2:blocksize(bind):blocksize(bind)^2 + (cind)*blocksize(bind)^2) = -data;
                                rind = rind + 1;                            
                            elseif(length(data) == blocksize(bind)^2)
                                %Create matrix and transpose it
                                temp = reshape(data,blocksize(bind),blocksize(bind))';
                                sdcone{blockindex(bind)}(:,cind+1) = -temp(:);
                                rind = rind + blocksize(bind);
                            else
                                error('Data read in dense entry A%d matrix is not the correct length. Expected first row of a matrix, or the entire matrix (all rows and columns, row major).',cind);
                            end
                        end
                    end
                    %Increment pointers
                    if(blocksize(bind) < 0 || rind > blocksize(bind))
                        rind = 1;
                        bind = bind + 1;
                    end
                    if(bind > N)
                        bind = 1;
                        cind = cind + 1;
                    end
                end
            case sINTEGER
                if(isempty(xtype))
                    error('Your integer declaration must come last!');
                end
                idx = sscanf(str,'*%d');
                xtype(idx) = 'I';
                
        end    
        eof = feof(fid);
    end
catch ME
    fclose(fid);
    rethrow(ME);
end

%Ensure sparse return args
if(~isempty(sdcone))
    for i = 1:length(sdcone)
        if(~issparse(sdcone{i}))
            sdcone{i} = sparse(sdcone{i});
        end
    end
end
if(~isempty(A))
    if(~issparse(A)), A = sparse(A); end    
end
if(print)
    nz = 0;
    if(~isempty(sdcone))
        for i = 1:length(sdcone)            
            nz = nz + nnz(sdcone{i});
        end
    end
    fprintf('Total Semidefinite NZs: %6d\n',nz);
    fprintf('Finished Reading SDPA File\n\n'); 
end
%Close
fclose(fid);

