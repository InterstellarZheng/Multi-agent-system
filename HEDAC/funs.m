function [value] = funs(x,z1,z2,z3,z4,N,T,delta_t,delta_x,m)%注意参数
    k = round(T/delta_t);
    Z1 = z1(1:(k+1),:);
    Z2 = z2(1:(k+1),:);
    Z3 = z3(1:(k+1),:);
    Z4 = z4(1:(k+1),:);
    value = subfuns(x,Z1,Z2,Z3,Z4,N,T,delta_t,delta_x,m)/intes(delta_x,Z1,Z2,Z3,Z4,N,T,delta_t,m);




function [value] = subfuns(x,z1,z2,z3,z4,N,T,delta_t,delta_x,m)
    compare = zeros(1,2);
    compare(1,1) = m-func(x,z1,z2,z3,z4,N,T,delta_t,delta_x);
    rhs = max(compare);
    value = rhs*rhs;



function [value] = intes(delta_x,z1,z2,z3,z4,N,T,delta_t,m)
    %针对不同区域需要更改
    xinit = [-0.5,-0.5];
    xfinal = [0.5,0.5];
    subvalue=0;
    while xinit(1,1) <= xfinal(1,1) %做积分运算
        while xinit(1,2) <= xfinal(1,2)
            x = xinit;
            subvalue = subvalue+subfuns(x,z1,z2,z3,z4,N,T,delta_t,delta_x,m);
            xinit(1,2)=xinit(1,2)+delta_x;
        end
        xinit(1,1)=xinit(1,1)+delta_x;
    end
    value = subvalue;
    
    
function [value] = func(x,z1,z2,z3,z4,N,T,delta_t,delta_x)%用于求解s
    value = subc(x,z1,z2,z3,z4,N,T,delta_t)/intec(delta_x,z1,z2,z3,z4,N,T,delta_t);
    

%用于求解c函数
function [value] = subc(x,z1,z2,z3,z4,N,T,delta_t)
    ktime = T/delta_t;
    value1=0;
    value2=0;
    value3=0;
    value4=0;
    
    for k = 1:ktime+1
        value1 = value1+RBF(x-z1(k,:))*delta_t;                 
        value2 = value2+RBF(x-z2(k,:))*delta_t;
        value3 = value3+RBF(x-z3(k,:))*delta_t;
        value4 = value4+RBF(x-z4(k,:))*delta_t;
    end
    if T == 0
        value = (value1+value2+value3+value4)/N;
    else
        value = (value1+value2+value3+value4)/(N*T);
    end
        
    
%用于求解c函数after integral
function [value] = intec(delta_x,z1,z2,z3,z4,N,T,delta_t)
    %针对不同区域需要更改
    xinit = [-0.5,-0.5];
    xfinal = [0.5,0.5];
    subvalue=0;
    while xinit(1,1) <= xfinal(1,1)
        while xinit(1,2) <= xfinal(1,2)
            x = xinit;
            subvalue = subvalue+subc(x,z1,z2,z3,z4,N,T,delta_t);
            xinit(1,2)=xinit(1,2)+delta_x;
        end
        xinit(1,1)=xinit(1,1)+delta_x;
    end
    value = subvalue;





    
    
    
    
    
    
    
    