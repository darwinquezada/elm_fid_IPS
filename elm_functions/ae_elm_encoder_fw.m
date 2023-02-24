function [B,ZH,H,RMSE,W,bmatrix] = ae_elm_encoder_fw(X,actfunc,numneurons,W,bias)
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
%   Hnew         -  Hidden layers

X = datanorm(X);    % Data normalization [0 1]
% % dlmwrite('originalData.csv',X,'delimiter',',','precision',7); 
% alph = size(X);
% % W=rand(numneurons,alph(2))*2-1;   % Input weights
%                                     
% % Orthogonal random weight
% wrow = numneurons;
% wcol = alph(2);
% wran = randn(wrow, wcol);
% W = orth( wran.' ).';               % Orthogonal random Matrix

bias = rand(numneurons,1);
bias = orth(bias.' ).';

Hprev = W*X' + bias;
ZH = actfunction(Hprev,actfunc);
B = pinv(ZH') * X;                  % Output weights Moore–Penrose
H=X*B';

% Get the RMSE
[RMSE, data_decoder] = ae_elm_decoder(X,H,B); % Decoder
end