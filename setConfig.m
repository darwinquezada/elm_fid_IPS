function conf = setConfig()
% Configuration

%% Configuration Data
conf.dataRepresentation = 'exponential';   % Data representation
                                        % - 'positive'
                                        % - 'normalized01'
                                        % - 'exponential'
                                        % - 'powed'
                                        % - 'raw'
conf.fullDataset        = './databases/MAN1.mat'; % Path .mat (Training and Test)

%% Data Positioning specifications
conf.distance           = 'meters';
conf.coordinates        = 'Local';       % Local coordinates 'Local'
                                         % WGS84
conf.distFunction       = 'cityblock';  % euclidean
                                         % cityblock
                                         % Manhattan
                                         % Minkowsky
                                         % Neyman
                                         % Soergel
                                         % cramariucPLGD40
                                         % cramariucPLGD10
                                         % squaredeuclidean

%% k-NN
conf.k                  = 11;

%% Auto-Encoder Extreme Learning Machine
conf.neurons            = 15;           % Number of Neurons > 0
conf.activationFunction = 'sig';        % Activation Function (We can Add more functions)
                                        %   'sig'       - Activation function sigmod
                                        %   'sin'       - Sine
                                        %   'tanh'      - Tangent Hyperbolic
end