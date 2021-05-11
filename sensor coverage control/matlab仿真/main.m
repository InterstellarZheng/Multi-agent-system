clear all;
close all;
delta_t=0.01;

T = 4000;
%定义初始变量
D = 5.0;
alpha1 = 0.01;
alpha2 = 0.001;
alpha3 = 32;

delta = 1;
V = D/delta;
alphak = -0.5;

%初始化
s1 = zeros(T+1,2);
s1(1,:) = [0,0];
s2 = zeros(T+1,2);
s2(1,:) = [0,1];
s3 = zeros(T+1,2);
s3(1,:) = [0,2];
s4 = zeros(T+1,2);
s4(1,:) = [0,3];
s5 = zeros(T+1,2);
s5(1,:) = [0,4];
s6 = zeros(T+1,2);
s6(1,:) = [0,5];
control1 = zeros(1,2);
control2 = zeros(1,2);
control3 = zeros(1,2);
control4 = zeros(1,2);
control5 = zeros(1,2);
control6 = zeros(1,2);





for k = 1:T
    
    agentposition = [s1(k,:);s2(k,:);s3(k,:);s4(k,:);s5(k,:);s6(k,:)];
    
%     [control1(:,1),control1(:,2)] = control_fun(s1(k,1),s1(k,2),s1(k,1),s1(k,2),s2(k,1),s2(k,2),s3(k,1),s3(k,2),s4(k,1),s4(k,2),s5(k,1),s5(k,2),s6(k,1),s6(k,2),delta,V,D);
%     [control2(:,1),control2(:,2)] = control_fun(s2(k,1),s2(k,2),s1(k,1),s1(k,2),s2(k,1),s2(k,2),s3(k,1),s3(k,2),s4(k,1),s4(k,2),s5(k,1),s5(k,2),s6(k,1),s6(k,2),delta,V,D);
%     [control3(:,1),control3(:,2)] = control_fun(s3(k,1),s3(k,2),s1(k,1),s1(k,2),s2(k,1),s2(k,2),s3(k,1),s3(k,2),s4(k,1),s4(k,2),s5(k,1),s5(k,2),s6(k,1),s6(k,2),delta,V,D);
%     [control4(:,1),control4(:,2)] = control_fun(s4(k,1),s4(k,2),s1(k,1),s1(k,2),s2(k,1),s2(k,2),s3(k,1),s3(k,2),s4(k,1),s4(k,2),s5(k,1),s5(k,2),s6(k,1),s6(k,2),delta,V,D);
%     [control5(:,1),control5(:,2)] = control_fun(s5(k,1),s5(k,2),s1(k,1),s1(k,2),s2(k,1),s2(k,2),s3(k,1),s3(k,2),s4(k,1),s4(k,2),s5(k,1),s5(k,2),s6(k,1),s6(k,2),delta,V,D);
%     [control6(:,1),control6(:,2)] = control_fun(s6(k,1),s6(k,2),s1(k,1),s1(k,2),s2(k,1),s2(k,2),s3(k,1),s3(k,2),s4(k,1),s4(k,2),s5(k,1),s5(k,2),s6(k,1),s6(k,2),delta,V,D);

    [control1(:,1),control1(:,2)] = control_fun(s1(k,1),s1(k,2),agentposition,delta,V,D);
    [control2(:,1),control2(:,2)] = control_fun(s2(k,1),s2(k,2),agentposition,delta,V,D);
    [control3(:,1),control3(:,2)] = control_fun(s3(k,1),s3(k,2),agentposition,delta,V,D);
    [control4(:,1),control4(:,2)] = control_fun(s4(k,1),s4(k,2),agentposition,delta,V,D);
    [control5(:,1),control5(:,2)] = control_fun(s5(k,1),s5(k,2),agentposition,delta,V,D);
    [control6(:,1),control6(:,2)] = control_fun(s6(k,1),s6(k,2),agentposition,delta,V,D);
    
    s1(k+1,:) = s1(k,:)+alphak.*control1(1,:);
    s2(k+1,:) = s2(k,:)+alphak.*control2(1,:);
    s3(k+1,:) = s3(k,:)+alphak.*control3(1,:);
    s4(k+1,:) = s4(k,:)+alphak.*control4(1,:);
    s5(k+1,:) = s5(k,:)+alphak.*control5(1,:);
    s6(k+1,:) = s6(k,:)+alphak.*control6(1,:);
    
    
end



%函数测试
% cob = zeros(1,2);
% [cob(:,1),cob(:,2)] = coortrans_gen_to_narr(1,2,5,4,0.1);


result = [s1(T+1,1),s1(T+1,2);
    s2(T+1,1),s2(T+1,2);
    s3(T+1,1),s3(T+1,2);
    s4(T+1,1),s4(T+1,2);
    s5(T+1,1),s5(T+1,2);
    s6(T+1,1),s6(T+1,2)];
figure(1);
scatter(result(:,1),result(:,2))
axis([-20,20,0,40]);
%title('Position')
%xlabel('xx');   ylabel('xy')
%hm = legend('无人机1','无人机2','无人机3','无人机4');
%set(hm,'Box','off');

figure(2);
plot(s1(:,1),s1(:,2))


















