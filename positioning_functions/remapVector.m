function [ m1 ] = remapVector( m0 , vin, vout )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
        
    m1 = m0*0;

    for i = 1:size(vin,2)
        m1(m0==vin(i))     = vout(i);
    end
    
end
