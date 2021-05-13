clear all;
close all;
%1x1正方形区域，中心为原点
delta_t=0.1;
delta_x = 0.05;
N = 4;
m = 1;%目标覆盖率，可更改

alpha = 0.002;
beita = 0.002;
gama = 0.5;
v=0.2;
%定义初始变量

%方便展示偏微分方程，对系数进行简化
k = beita/alpha;
n = delta_x*delta_x*k+4;

%求解偏微分方程
xnumber = 1/delta_x;    %1对应1个单位长度，一段长度需要两个数据点
%tnumber = 20/delta_t;   %20对应20s，修改20即为修改运行时间
tnumber = 100;
A = zeros(xnumber+1,xnumber+1,tnumber);  %+1是为了将初始条件也放到矩阵当中



s1 = zeros(tnumber,2);%初始位置t0为矩阵第一行
s1(1,:) = [0.25,0.25];
s2 = zeros(tnumber,2);
s2(1,:) = [-0.25,0.25];
s3 = zeros(tnumber,2);
s3(1,:) = [-0.25,-0.25];
s4 = zeros(tnumber,2);
s4(1,:) = [0.25,-0.25];




for i = 2:xnumber
    for j = 2:xnumber
        A(i,j,1) = 0.01*exp((-1)*(((i-1)*delta_x-0.5)*((i-1)*delta_x-0.5)+((j-1)*delta_x-0.5)*((j-1)*delta_x-0.5)));             %-0.5为了将矩阵坐标系转换为实际的坐标系
    end
end
%边界法向量为0，相邻两个数据相等
A(1,2:xnumber,1) = A(2,2:xnumber,1);
A(xnumber+1,2:xnumber,1) = A(xnumber,2:xnumber,1);
A(1:xnumber+1,1,1) = A(1:xnumber+1,2,1);
A(1:xnumber+1,xnumber+1,1) = A(1:xnumber+1,xnumber,1);



    




for t = 2:tnumber           %t为下一时间点
    T = (t-2)*delta_t;  %t=1时实际T为0
    
    for i = 2:xnumber
        for j = 3:xnumber
            x = [(i-1)*delta_x-0.5,(j-1)*delta_x-0.5];%转换为实际的平面xy坐标
            f = (delta_x*delta_x)*((gama/alpha)*funa(x,s1,s2,s3,s4,T,delta_t,delta_x)-funs(x,s1,s2,s3,s4,N,T,delta_t,delta_x,m));
            A(i,j,t) = n*A(i,j-1,t-1)-A(i+1,j-1,t-1)-A(i-1,j-1,t-1)-A(i,j-2,t-1)+f;                       %t可能为0因此考虑是否存在时间导数
            if i == 2
                A(1,j,t) = A(2,j,t);
            end
            if i == xnumber
                A(xnumber+1,j,t) = A(xnumber,j,t);
            end
        end
    end
    A(1:xnumber+1,xnumber+1,t) = A(1:xnumber+1,xnumber,t);
        
        dudx1 = (A(floor((s1(t-1,1)+0.5)/delta_x)+1,floor((s1(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s1(t-1,1)+0.5)/delta_x)+2,floor((s1(t-1,2)+0.5)/delta_x)+1,t))/delta_x;
        dudy1 = (A(floor((s1(t-1,1)+0.5)/delta_x)+1,floor((s1(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s1(t-1,1)+0.5)/delta_x)+1,floor((s1(t-1,2)+0.5)/delta_x)+2,t))/delta_x;
        
        dudx2 = (A(floor((s2(t-1,1)+0.5)/delta_x)+1,floor((s2(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s2(t-1,1)+0.5)/delta_x)+2,floor((s2(t-1,2)+0.5)/delta_x)+1,t))/delta_x;
        dudy2 = (A(floor((s2(t-1,1)+0.5)/delta_x)+1,floor((s2(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s2(t-1,1)+0.5)/delta_x)+1,floor((s2(t-1,2)+0.5)/delta_x)+2,t))/delta_x;
        
        dudx3 = (A(floor((s3(t-1,1)+0.5)/delta_x)+1,floor((s3(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s3(t-1,1)+0.5)/delta_x)+2,floor((s3(t-1,2)+0.5)/delta_x)+1,t))/delta_x;
        dudy3 = (A(floor((s3(t-1,1)+0.5)/delta_x)+1,floor((s3(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s3(t-1,1)+0.5)/delta_x)+1,floor((s3(t-1,2)+0.5)/delta_x)+2,t))/delta_x;
        
        dudx4 = (A(floor((s4(t-1,1)+0.5)/delta_x)+1,floor((s4(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s4(t-1,1)+0.5)/delta_x)+2,floor((s4(t-1,2)+0.5)/delta_x)+1,t))/delta_x;
        dudy4 = (A(floor((s4(t-1,1)+0.5)/delta_x)+1,floor((s4(t-1,2)+0.5)/delta_x)+1,t)...
            -A(floor((s4(t-1,1)+0.5)/delta_x)+1,floor((s4(t-1,2)+0.5)/delta_x)+2,t))/delta_x;
        
        dudx1 = dudx1/sqrt(dudx1*dudx1+dudy1*dudy1);
        dudy1 = dudy1/sqrt(dudx1*dudx1+dudy1*dudy1);
        dudx2 = dudx2/sqrt(dudx2*dudx2+dudy2*dudy2);
        dudy2 = dudy2/sqrt(dudx2*dudx2+dudy2*dudy2);
        dudx3 = dudx3/sqrt(dudx3*dudx3+dudy3*dudy3);
        dudy3 = dudy3/sqrt(dudx3*dudx3+dudy3*dudy3);
        dudx4 = dudx4/sqrt(dudx4*dudx4+dudy4*dudy4);
        dudy4 = dudy4/sqrt(dudx4*dudx4+dudy4*dudy4);
        disp(dudx1*v*delta_t);
        s1(t,1) = s1(t-1,1)+dudx1*v*delta_t;
        s1(t,2) = s1(t-1,2)+dudy1*v*delta_t;
        s2(t,1) = s2(t-1,1)+dudx2*v*delta_t;
        s2(t,2) = s2(t-1,2)+dudy2*v*delta_t;
        s3(t,1) = s3(t-1,1)+dudx3*v*delta_t;
        s3(t,2) = s3(t-1,2)+dudy3*v*delta_t;
        s4(t,1) = s4(t-1,1)+dudx4*v*delta_t;
        s4(t,2) = s4(t-1,2)+dudy4*v*delta_t;
%         if s1(t,1) > 0.5
%             s1(t,1) = 0.5;
%         end
%         if s1(t,1) < -0.5
%             s1(t,1) = -0.5;
%         end
%         if s1(t,2) > 0.5
%             s1(t,2) = 0.5;
%         end
%         if s1(t,2) < -0.5
%             s1(t,2) = -0.5;
%         end
%         
%         if s2(t,1) > 0.5
%             s2(t,1) = 0.5;
%         end
%         if s2(t,1) < -0.5
%             s2(t,1) = -0.5;
%         end
%         if s2(t,2) > 0.5
%             s2(t,2) = 0.5;
%         end
%         if s2(t,2) < -0.5
%             s2(t,2) = -0.5;
%         end
%         
%         if s3(t,1) > 0.5
%             s3(t,1) = 0.5;
%         end
%         if s3(t,1) < -0.5
%             s3(t,1) = -0.5;
%         end
%         if s3(t,2) > 0.5
%             s3(t,2) = 0.5;
%         end
%         if s3(t,2) < -0.5
%             s3(t,2) = -0.5;
%         end
%         
%         if s4(t,1) > 0.5
%             s4(t,1) = 0.5;
%         end
%         if s4(t,1) < -0.5
%             s4(t,1) = -0.5;
%         end
%         if s4(t,2) > 0.5
%             s4(t,2) = 0.5;
%         end
%         if s4(t,2) < -0.5
%             s4(t,2) = -0.5;
%         end
        disp(t);
end

plot(s1(:,1),s1(:,2),s2(:,1),s2(:,2),s3(:,1),s3(:,2),s4(:,1),s4(:,2))



