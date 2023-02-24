function out = actfunction(data, type)
% [out] = activation_function(data, type)
% Activation function 
%   data    :   Vector or matrix 
%   type    :
%   'sig'       - Activation function sigmod
%   'sin'       - Sine
%   'tanh'      - Tangent Hyperbolic
% Out       :   Activation function output

    switch lower(type)
        case {'relu'}
            out = relu(data);
        case {'elu'}
            out = elulayer(data);
        case {'leakyRelu'}
            out = leakyrelu(data);
        case {'sig','sigmoid'}
            out = 1 ./ (1 + exp(-data));
        case {'sin','sine'}
            out = sin(data);
        case {'cos','cosine'}
            out = cos(data);
        case {'tanh'}
            out = (exp(data)-exp(-data))/(exp(data)+exp(-data)); 
        case {'tansig'}
            out = 2./(1+exp(-2*data))-1;
    end
end