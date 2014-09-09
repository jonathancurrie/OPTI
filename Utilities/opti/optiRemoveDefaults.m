function opts = optiRemoveDefaults(opts,defs)
%OPTIREMOVEDEFAULTS  Remove Default Options from an Options Structure
try
    oFn = fieldnames(opts);
    for i = 1:length(oFn)
        label = oFn{i};
        if(isfield(defs,label))
            if(ischar(opts.(label)))
                if(strcmpi(defs.(label),opts.(label)))
                    opts = rmfield(opts,label);
                end
            else
                if(defs.(label) == opts.(label))
                    opts = rmfield(opts,label);
                end
            end
        end
    end
catch ME %do nothing, accept defaults [J.Lofberg]
    optiwarn('OPTI:DefaultOpts',['There was an error removing the default options, you may experience significantly slower solver performance. '...
                                 'This is typically due to an old version of MATLAB (pre 2010a). The actual error is listed below:\n\n%s\n'],ME.message);
end