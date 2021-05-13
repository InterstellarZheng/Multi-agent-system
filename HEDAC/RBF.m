%用于求解RBF函数值
function [value] = RBF(deltax)
    value = exp((-1)*(deltax(1,1)*deltax(1,1)+deltax(1,2)*deltax(1,2))/2)/(2*pi);
    

