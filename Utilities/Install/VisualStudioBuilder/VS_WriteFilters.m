function VS_WriteFilters(projPath,projName,upath,source,header)

global fh;

%Header
docNode = com.mathworks.xml.XMLUtils.createDocument('Project');
p = docNode.getDocumentElement;
p.setAttribute('ToolsVersion','4.0');
p.setAttribute('xmlns','http://schemas.microsoft.com/developer/msbuild/2003');

%Generate GUIDs
cdir = cd;
cd(GETCD);
dos(['genGUID ' num2str(size(source,1) + size(header,1) + 3)]);
cd(cdir);
%Open Text File
fh = fopen([GETCD 'guids.txt']);

%Base Filter Identifiers
pc = docNode.createElement('ItemGroup');
pc.appendChild(createSource(docNode,[]));
pc.appendChild(createHeader(docNode,[]));
pc.appendChild(createResource(docNode,[]));

%Source Filters
p.appendChild(createCustomFilters(docNode,pc,upath,source,'source'));
%Header Filters
p.appendChild(createCustomFilters(docNode,pc,upath,header,'header'));

%Add Source Files to Filters
p.appendChild(createFileList(docNode,upath,source,'ClCompile','Source Files'));
%Add Header Files to Filters
p.appendChild(createFileList(docNode,upath,header,'ClInclude','Header Files'));

%Clean Up Empty Filters
VS_CleanFilters(docNode);

%Write Whole Project
xmlwrite([projPath '\' projName '_filters.xml'],docNode);
xmlwrite([projPath '\' projName '.vcxproj.filters'],docNode);

%Close GUID folder
fclose(fh);


function pc = createFileList(docNode,upath,files,elem,filter)
pc = docNode.createElement('ItemGroup');
for i = 1:size(files,1) %for each path    
    if(iscell(upath))
        p = upath{files{i,3}};
    else
        p = upath;
    end
   
    len = length(p);
    if(any(strfind(p,':')))
        ind = strfind(p,'\');
        len2 = ind(end)+1;
    else
        len2 = 1;
    end    
    attp = files{i,1}(len2:end);
    txtp = files{i,1}(len+1:end);
    
    if(~isempty(files{i,2}))
        for j = 1:size(files{i,2},1) %for each file
            pc1 = docNode.createElement(elem);
            pc1.setAttribute('Include',['..\' attp '\' files{i,2}{j}]);
            addElemText(docNode,pc1,'Filter',[filter txtp]);
            pc.appendChild(pc1);
        end
    end
end

function pc = createCustomFilters(docNode,pc,upath,files,mode)
for i = 1:size(files,1)
    s = files{i,1};
    if(~isempty(s))    
        if(iscell(upath))
            p = upath{files{i,3}};
            len = length(p);
            ind = strfind(s,p);
        else
            len = length(upath);
            ind = strfind(s,upath);
        end
        s(ind:len(1)) = [];
        switch(mode)
            case 'source'
                pc.appendChild(createSource(docNode,s));
            case 'header'
                pc.appendChild(createHeader(docNode,s));
            otherwise
                error('unknown filter mode');
        end
    end
end
        
function pc = createSource(docNode,path)
pc = docNode.createElement('Filter');
if(~isempty(path))
    pc.setAttribute('Include',['Source Files' path]);
else
    pc.setAttribute('Include','Source Files');
end
addElemText(docNode,pc,'UniqueIdentifier',randIdent());
addElemText(docNode,pc,'Extensions','cpp;c;cc;cxx;def;odl;idl;hpj;bat;asm;asmx');

function pc = createHeader(docNode,path)
pc = docNode.createElement('Filter');
if(~isempty(path))
    pc.setAttribute('Include',['Header Files' path]);
else
    pc.setAttribute('Include','Header Files');
end
addElemText(docNode,pc,'UniqueIdentifier',randIdent());
addElemText(docNode,pc,'Extensions','h;hpp;hxx;hm;inl;inc;xsd');

function pc = createResource(docNode,path)
pc = docNode.createElement('Filter');
if(~isempty(path))
    pc.setAttribute('Include',['Resource Files' path]);
else
    pc.setAttribute('Include','Resource Files');
end
addElemText(docNode,pc,'UniqueIdentifier',randIdent());
addElemText(docNode,pc,'Extensions','rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe;resx;tiff;tif;png;wav;mfcribbon-ms');


function addElemText(docNode,pc,elem,text)
c = docNode.createElement(elem);
c.appendChild(docNode.createTextNode(text));
pc.appendChild(c);


function str = randIdent() %reads from pre-generated file

global fh;
str = fgetl(fh);
if(length(str) < 5)
    error('Error reading GUID');
end

%Get cwd of this file
function str = GETCD()
str = which('VS_WriteProj.m');
ind = strfind(str,'\');
str(ind(end)+1:end) = [];