function VS_CleanFilters(docNode)

doc     = docNode.getChildNodes;    
proj    = doc.item(0).getChildNodes; 
filt    = proj.item(0).getChildNodes; %filter list
src     = proj.item(1).getChildNodes; %source file list
hdr     = proj.item(2).getChildNodes; %header file list

%Build Complete File List Saving only Filter Locations
len = src.getLength + hdr.getLength;
ind = 1;
files = cell(len,1);
%Read in source files
for i = 1:src.getLength
    srcC = src.item(i-1).getChildNodes;
    files{ind,1} = char(srcC.getChildNodes.getTextContent);
    ind = ind + 1;
end
%Read in header files
for i = 1:hdr.getLength
    hdrC = hdr.item(i-1).getChildNodes;
    files{ind,1} = char(hdrC.getChildNodes.getTextContent);
    ind = ind + 1;
end

%For each Filter added check if we have a file present or subdirectory
%present, if not, delete the child node
for i = 4:filt.getLength %skip default first 3
    if(i > filt.getLength)
        break;
    end
    f = filt.item(i-1).getChildNodes;
    att = f.getAttributes;
    p = char(att.item(0).getValue);
    if(all(cellfun(@isempty,strfind(files,p))))
        filt.removeChild(f);
    end
end