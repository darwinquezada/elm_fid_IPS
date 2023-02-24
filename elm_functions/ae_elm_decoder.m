function [RMSE, data_decoder] = ae_elm_decoder(X,H,B)
    data_decoder = H*pinv(B');
    RMSE = sqrt(mse(X,data_decoder)); %Problem here -- let's check later
end