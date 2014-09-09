function [str,deflist] = readMissingSymbols(file)
%Reads a file like:
%lpsolve.obj : error LNK2019: unresolved external symbol _create_hash_table referenced in function _mainloop 
%lpsolve.obj : error LNK2019: unresolved external symbol _free_hash_table referenced in function _ExitFcn 
%lpsolve.obj : error LNK2019: unresolved external symbol _findhash referenced in function _impl_get_handle 
%
% And returns a guess at the missing symbol names

data = fileread(file);

st = 'external symbol';
en = 'referenced in';
idx1 = strfind(data,st) + length(st);
idx2 = strfind(data,en) - 2;

if(length(idx1) ~= length(idx2))
    error('unequal start and end vectors');
end

str = cell(length(idx1),1);

for i = 1:length(idx1)
    str{i} = processSymbol(data(idx1(i):idx2(i)));
end

%Remove copies 
% str = unique(str);

%Build a list for definition file
deflist = '[';
for i = 1:length(str)-1
    deflist = sprintf('%s''%s'',',deflist,str{i});
end
deflist = sprintf('%s''%s'']',deflist,str{end});


function str = processSymbol(symb)
%Find start
for i = 1:length(symb)
    if(~isspace(symb(i)) && symb(i) ~= '_')
        break;
    end
end
%Find end 
str = symb(i:end);
idx = strfind(str,'@');
if(~isempty(idx))
    str = str(1:idx(1)-1);
end




    

