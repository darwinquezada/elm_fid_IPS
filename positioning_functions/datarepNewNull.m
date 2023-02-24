function [ m1 ] = datarepNewNull( m0 , oldNull, newNull )
    m1= m0.*(m0~=oldNull)+newNull.*(m0==oldNull);
end
