function data = dataRead(path)
    slitPath = split(path,'.');
    sz = size(slitPath,1);
    ext = slitPath{sz};
    % newdata = [];
    switch lower(ext)
        case 'mat'
           data = load(path);
        case 'csv'
           file = csvread(path,1,0);
           data       = [file];
        otherwise
            msg = 'File error. Files supported .mat and .csv';
            error(msg);
    end
    
    
end