function parameter = initializeUniform(sz,bound)

Z = 2*rand(sz,'single') - 1;
parameter = bound * Z;
% parameter = dlarray(parameter);

end