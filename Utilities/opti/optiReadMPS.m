function prob = optiReadMPS(filename)
%OPTIREADMPS  Read an MPS/QPS file and convert to Matlab Values
%
% prob = optiReadMPS(filename) reads the file specified by filename and
% converts it to an optiprob. If you do not specify a file a extension, it
% will default to .mps, otherwise you can read .mps or .qps files using
% this interface. The returned structure is solver independent so the user
% can manually extract matrices as required.
%
% MPS & Parser Assumptions made:
%   - ALL names are CASE SENSITIVE
%   - Default bounds are 0 <= x <= inf
%   - Integer constraints are NOT binary. To set a binary constraint, use
%   the bound BV or set an upper bound of 1 (i.e. using UI)
%   - Comments (*) are only at the top of the file
%   - Maximum of 7 columns (fields) for section COLUMNS + RHS
%   - QUADOBJ only uses 3 columns (i, j, val)
%   - SOS & RANGES not currently implemented
%
% Note this is particularly inefficient code written quickly - you are best
% to use coinRead() instead! This function will be removed in a future
% release.

%   Copyright (C) 2011 Jonathan Currie (Control Engineering)

% Make sure we have file extension!
if(isempty(strfind(filename,'.mps')) && isempty(strfind(filename,'.qps')) && isempty(strfind(filename,'.')))
    filename = [filename '.mps'];
end

%Open File
fid = fopen(filename);
if(fid == -1)
    error('Error opening file %s',filename);
end

%Default Internals
sense = -1;
H = []; lb = []; ub = []; ivarsB = [];
eof = 0;

try
    %Main Loop
    while(~eof)
        %Read MPS File String
        str = fgetl(fid);    
        if(any(str) && ~isblank(str) && ~isComment(str)) %check for header           
            %Act depending on heading
            str = textscan(str,'%s %s'); %extract just the header
            heading = str{1}{1};
            switch(lower(heading))
                case 'name'
                    if(isempty(str{2}))
                        name = '';
                    else
                        name = str{2}{1};
                    end                
                case 'objsense'
                    str = fscanf(fid,'%s',1);
                    switch(lower(str))
                        case 'max'
                            sense = 1;
                        case 'min'
                            sense = -1;
                        otherwise
                            error('Unknown objsense: %s',str);
                    end   
                case 'objname'
                    continue; %ignore
                case 'rows'
                    [CON,e,COST] = readROWS(fid);
                case 'columns' 
                    [f,A,DECVAR,ivars] = readCOLUMNS(fid,CON,COST);
                case 'rhs'
                    b = readRHS(fid,CON);
                case 'ranges'
                    error('RANGES not implemented yet');
                case 'bounds'
                    [lb,ub,ivarsB] = readBOUNDS(fid,DECVAR);
                case 'quadobj'
                    H = readQUADOBJ(fid,DECVAR);
                case 'sos'
                    error('SOS not implemented yet');
                case 'endata'
                    break;
                otherwise
                    error('Unknown MPS header: %s',heading);
            end
        end
        eof = feof(fid);
    end
catch ME
    fclose(fid);
    rethrow(ME);
end
    
%Close
fclose(fid);
%Check we have bounds
if(isempty(lb)), lb = zeros(length(DECVAR),1); end
if(isempty(ub)), ub = inf*ones(length(DECVAR),1); end
%Correct Sense
f = f*-sense;
%Work out Integer Vars
if(~isempty(ivarsB))
    ivar = find((ivars + ivarsB) > 0);
else
    ivar = find(ivars > 0);
end

%Build return object
prob = optiprob('name',name,'H',H,'f',f,'mix',A,b,e,'bounds',lb,ub,'int',ivar);


function [CON,e,COST] = readROWS(fid)
%Read the section ROWS creating a container of constraints and cost
%function

CON = containers.Map(); nCon = int32(1); 
COST = containers.Map(); nCost = int32(1);
fpos = ftell(fid); %save file position for coming back here after determing # rows
eof = 0; 
nRows = 0;

%Determine # lines to read
while(~eof)
    str = fgetl(fid);    
    if(any(str) && ~isblank(str)) %header
        break;
    end    
    if(str)
        nRows = nRows + 1;
    end
    %Check eof
    eof = feof(fid);
end
%return to ROWS
fseek(fid,fpos,-1);
%Read Rows into cell array
C = textscan(fid,'%s %s',nRows);
%Preallocate
e = zeros(nRows-1,1); %ASSUME only one N
%Process Rows for Constraints + Cost
for i = 1:nRows    
    str = C{1}{i};
    var = C{2}{i};
    if(~isempty(str))
        switch(lower(str))
            case 'l'
                CON(var) = nCon; 
                e(nCon) = -1;
                nCon = nCon + 1;
            case 'g'
                CON(var) = nCon; 
                e(nCon) = 1;
                nCon = nCon + 1;
            case 'e'
                CON(var) = nCon; 
                e(nCon) = 0;
                nCon = nCon + 1;
            case 'n'
                COST(var) = nCost; 
                nCost = nCost + 1;
        end
    end
end
%Check Sizes
if(nCost > 2) %+1 above, so should be 1
    error('Expected only one objective function! - Found %d N entries in ROWS',nCost-1);
end
                

function [f,A,DECVAR,ivars] = readCOLUMNS(fid,CON,COST)
%Read the section COLUMNS returning the f and A matrices + integer vars

DECVAR = containers.Map; nDec = int32(1);
MARKER = containers.Map; nMar = int32(1);

fpos = ftell(fid); %save file position for coming back here after determing lengths
eof = 0; 
nzA = 0;
intOn = 0;
nRows = 0;

%Determine # lines to read
while(~eof)
    str = fgetl(fid);    
    if(any(str) && ~isblank(str)) %header
        break;
    end    
    if(str)
        nRows = nRows + 1;
    end
    %Check eof
    eof = feof(fid);
end
%return to COLUMNS
fseek(fid,fpos,-1);
%Read Columns into cell array
C = textscan(fid,'%s %s %s %s %s %s %s',nRows);
%Process Columns for Sizes + Decision Vars
for i = 1:nRows
    for j = [1 2 4 6]
        str = C{j}{i};
        if(~isempty(str))
            switch(j)
                case 1 %dec var or marker                    
                    if(strcmp(C{2}{i},'''MARKER'''))
                        MARKER(str) = nMar; nMar = nMar + 1;
                    elseif(~isKey(DECVAR,str))
                        DECVAR(str) = nDec; nDec = nDec + 1;
                    end
                case {2,4,6} %con var
                   if(isKey(CON,str))
                       nzA = nzA + 1;
                   end
            end
        end
    end    
end
%Preallocate sparse triples + cost
Ai = zeros(nzA,1); ai = 1;
Aj = zeros(nzA,1);
As = zeros(nzA,1);
f = zeros(length(DECVAR),1); skip = 0;
ivars = zeros(length(DECVAR),1);
%Read in Data
for i = 1:nRows
    skipj = 0;
    for j = 1:7    
        str = C{j}{i};
        if(~isempty(str))
            switch(j-skipj)
                case 1 %dec var or Marker
                    if(isKey(DECVAR,str))
                        m = DECVAR(str);
                        if(intOn)
                            ivars(m) = 1;
                        end
                    elseif(isKey(MARKER,str))
                        if(strcmp(C{3}{i},'''INTORG'''))
                            intOn = 1;
                        elseif(strcmp(C{3}{i},'''INTEND'''))
                            intOn = 0;
                        else
                            error('Unknown Marker Argument: %s',C{3}{i});
                        end
                        skipj = 7; %skip remainder of j indices
                    else
                        error('Unknown argument (expected dec var or marker): %s',str);
                    end                    
                case {2,4,6} %con or cost var
                    if(isKey(CON,str))
                        n = CON(str);
                    elseif(isKey(COST,str))
                        val = str2double(C{j+1}{i});
                        if(isnan(val)), error('Error Converting Cost Value %s',str); end
                        f(m) = val;
                        skip = 1; %skip next val (already read it)
                    else
                        error('Unknown argument: %s',str);
                    end
                case {3,5,7} %val
                    if(~skip)
                        val = str2double(C{j}{i});
                        if(isnan(val)), error('Error Converting Cons Value %s',str); end
                        Ai(ai) = m;
                        Aj(ai) = n;
                        As(ai) = val;
                        ai = ai + 1;
                    else
                        skip = 0;
                    end
            end
        end
    end    
end

%Check if we have 0 in index
if(any(Aj == 0))
    error('Failure reading columns (nnz != indicies saved) - check for errors!');
end

%Build A
A = sparse(Aj,Ai,As,length(CON),length(DECVAR)); %note transposed i,j


function b = readRHS(fid,CON)
%Read the section RHS return the b vector

eof = 0; 
nRows = 0;
skip = 0;
fpos = ftell(fid);
%Preallocate
b = zeros(length(CON),1);
%Determine # lines to read
while(~eof)
    str = fgetl(fid);    
    if(any(str) && ~isblank(str)) %header
        break;
    end    
    if(str)
        nRows = nRows + 1;
    end
    %Check eof
    eof = feof(fid);
end
%return to RHS
fseek(fid,fpos,-1);
%Read Rows into cell array
C = textscan(fid,'%s %s %s %s %s %s %s',nRows);
%Process Rows for Bounds
for i = 1:nRows  
    for j = 2:7 %ignore rhs name (assume only one!)
        str = C{j}{i};
        if(~isempty(str))
            switch(j)
                case {2,4,6} %cons
                    if(isKey(CON,str))
                        ind = CON(str);
                    else
                        optiwarn('opti:mps','Unknown Constraint when reading RHS: %s, skipping',str);
                        skip = 1;
                    end
                case {3,5,7} %val
                    if(~skip)
                        val = str2double(str);
                        if(isnan(val)), error('Cannot convert RHS value: %s',str); end
                        b(ind) = val;
                    else
                        skip = 0;
                    end
            end
        end
    end    
end

    


function [lb,ub,ivars] = readBOUNDS(fid,DECVAR)
%Read the section BOUNDS and return lb, ub + any other integer vars

fpos = ftell(fid); %save file position for coming back here after determing lengths

%Preallocate
lb = zeros(length(DECVAR),1); %MPS assumes 0 <= x < inf ??
ub = inf*ones(length(DECVAR),1);
ivars = zeros(length(DECVAR),1);
eof = 0; 
nRows = 0;

%Determine # lines to read
while(~eof)
    str = fgetl(fid);    
    if(any(str) && ~isblank(str)) %header
        break;
    end    
    if(str)
        nRows = nRows + 1;
    end
    %Check eof
    eof = feof(fid);
end
%return to BOUNDS
fseek(fid,fpos,-1);
%Read Rows into cell array
C = textscan(fid,'%s %s %s %s',nRows);
%Process Rows for Bounds
for i = 1:nRows    
    boundT = C{1}{i};
    decVar = C{3}{i};
    val = str2double(C{4}{i});
    if(isnan(val)), error('Could not convert Bound %s %s %s',boundT,decVar,val); end
    if(isKey(DECVAR,decVar))
        ind = DECVAR(decVar);
    else
        error('Unknown decision variable: %s',decVar);
    end
    if(~isempty(str))
        switch(lower(boundT))
            case 'up'
                ub(ind) = val;
            case 'lo'
                lb(ind) = val;
            case 'ui'
                ub(ind) = val;
                ivars(ind) = 1;
            case 'li'
                lb(ind) = val;
                ivars(ind) = 1;
            case 'fx'
                lb(ind) = val;
                ub(ind) = val;                    
            case 'bv'
                lb(ind) = 0;
                ub(ind) = 1;
                ivars(ind) = 1;
            case 'mi'
                lb(ind) = -inf;
            case 'pl'
                ub(ind) = inf;
            case 'fr'
                lb(ind) = -inf;
                ub(ind) = inf;
            otherwise
                error('Unknown Bound Type: %s',boundT);
        end
    end
end


function H = readQUADOBJ(fid,DECVAR)
%Read the section QUADOBJ and return H
%NOTE - I read quadobj has the same format as columns - but only ever seen
%3 columns in quadobj. Assuming for now this is the case!

fpos = ftell(fid); %save file position for coming back here after determing # rows
eof = 0; 
nRows = 0;

%Determine # lines to read
while(~eof)
    str = fgetl(fid);    
    if(any(str) && ~isblank(str)) %header
        break;
    end    
    if(str)
        nRows = nRows + 1;
    end
    %Check eof
    eof = feof(fid);
end
%Preallocate                
Qi = zeros(nRows,1); qi = 1;
Qj = zeros(nRows,1);
Qs = zeros(nRows,1);
%return to ROWS
fseek(fid,fpos,-1);
%Read Rows into cell array
C = textscan(fid,'%s %s %s',nRows);
%Process Rows for QuadObj
for i = 1:nRows    
    x1 = C{1}{i};
    x2 = C{2}{i};
    str = C{3}{i};
    if(~isempty(str))
        if(isKey(DECVAR,x1))
            m = DECVAR(x1);
        else
            error('Unknown Decision Variable when Reading QUADOBJ: %s',x1);
        end
        if(isKey(DECVAR,x2))
            n = DECVAR(x2);
        else
            error('Unknown Decision Variable when Reading QUADOBJ: %s',x1);
        end
        val = str2double(str);
        if(isnan(val)), error('Error converting quadobj val: %s',str); end
        Qi(qi) = m;
        Qj(qi) = n;
        Qs(qi) = val;
        qi = qi + 1;
    end
end

%Check if we have 0 in index
if(any(Qj == 0))
    error('Failure reading quadobj (nnz != indicies saved) - check for errors!');
end

%Build H
H = sparse(Qj,Qi,Qs,length(DECVAR),length(DECVAR)); %note transposed i,j


%Check if first character is blank (not a header)
function x = isblank(str)
if(str(1) == ' ')
    x = 1;
else
    x = 0;
end
    
%Check if first character is a star (comment)
function x = isComment(str)
if(str(1) == '*')
    x = 1;
else
    x = 0;
end   
 
