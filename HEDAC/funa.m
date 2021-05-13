function [value] = funa(x,z1,z2,z3,z4,T,delta_t,delta_x)
    k = round(T/delta_t);
    Z1 = z1(1:(k+1),:);
    Z2 = z2(1:(k+1),:);
    Z3 = z3(1:(k+1),:);
    Z4 = z4(1:(k+1),:);
    value = suba(x,Z1,Z2,Z3,Z4)/intea(delta_x,Z1,Z2,Z3,Z4);




function [value] = suba(x,z1,z2,z3,z4)
    value1=RBF(x-z1);
    value2=RBF(x-z2);
    value3=RBF(x-z3);
    value4=RBF(x-z4);
    value = value1+value2+value3+value4;



function [value] = intea(delta_x,z1,z2,z3,z4)
    %针对不同区域需要更改
    xinit = [-0.5,-0.5];
    xfinal = [0.5,0.5];
    subvalue=0;
    while xinit(1,1) <= xfinal(1,1)
        while xinit(1,2) <= xfinal(1,2)
            x = xinit;
            subvalue = subvalue+suba(x,z1,z2,z3,z4);
            xinit(1,2)=xinit(1,2)+delta_x;
        end
        xinit(1,1)=xinit(1,1)+delta_x;
    end
    value = subvalue;