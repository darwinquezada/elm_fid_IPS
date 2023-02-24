function [B,ZH,H,RMSE,W,bias] = ae_elm_encoder(X,actfunc,numneurons)
% [B,Hold,Hnew] = ae_elm_encoder(data,actfunction,numneurons)
% Autoencoder function using Extreme Learning Machine (ELM)
% Input parameters:
%   X            -  Data to encoder Matrix(m x n) or vector
%   actfunc      -  Activation function (help actfunct)
%                       'sig'
%                       'sin'
%                       'tanh'
%   numneurons   -  Number of neurons 
% Output function:
%   B            -  Beta B = [B1, B2, ..., Bn]' is the output weight matrix
%   ZH           -  Activation function result (compression)
%   H            -  Hidden layers

X = datanorm(X);    % Data normalization [0 1]
% dlmwrite('originalData.csv',X,'delimiter',',','precision',7); 
alph = size(X);
W=rand(numneurons,alph(2))*2-1;   % Input weights
W = orth( W.' ).';
                 
% Orthogonal random weight
%wrow = numneurons;
%wcol = alph(2);
%wran = randn(wrow, wcol);
%Wor = orth( wran.' ).';
%W = bsxfun(@rdivide, Wor, sqrt(sum(Wor.^2, 2))); % Orthogonal random Weight Matrix Thanks Matlab!!!

bias = rand(numneurons,1);
% bias = orth(bias.' ).';
bias = orth(bias);

Hprev = W*X' + bias;
ZH = actfunction(Hprev,actfunc);
B = pinv(ZH') * X;                  % Output weights Moore–Penrose
H=X*B';

% Get the RMSE
data_decoder = H*pinv(B');
RMSE = sqrt(mse(X,data_decoder));
%[RMSE, data_decoder] = ae_elm_decoder(X,H,B); % Decoder

end