function db = data_representation(X,type)
% data = dataRepresentation(X,type)
% X     : RSS
% Type  :   'positive', 
%           'exponential'
%           'powed'

    type = type(find(~isspace(type)));
    % minval = min(min(X(:)));
    % maxval = max(max(X(:)));
    
    switch lower(type)
        case {'positive'}
            % data = X + 101;
            % data(data == 201) = 0;
            db = datarepPositive( X );
        case 'exponential'
            db = datarepPositive( X );
        case {'powed'}
            db = datarepExponential( X );
        otherwise
            msg = 'Error: Incorrect type';
            error(msg);
    end
end