%clear;
%clc;
%function [F_simple,u_center] = calc_dragforce(Rp,tau,theta,n)
function [F_drag0,F_lift,u_uniform,CD] = calc_dragforce(Rp,tau,theta,n)

    %基础参数
    g = 9.8;
    Dp = 2*Rp;
    rhof = 1000;
    nv = 1e-6;
    shields = tau/1650/g/Dp;
    
    %升力系数直接定为定值
    CL = 0.178;

    %计算速度分布
    delta = 1.0*Dp;  %线性层厚度
    u_star = sqrt(tau/rhof);  %剪切流速
    ks = 2*Dp;
    ks_plus = u_star*ks/nv;  %粗糙程度
    C = 0;
    if ks_plus <= 1000
        C = -0.993*log(ks_plus) + 12.36;
    else
        C = 5.5;
    end
    
    %速度波动的标准差
    sigma_u = 2*u_star;

    %速度分布函数
    u_y = @(y) u_star*C/delta.*y.*(y <= delta) + ...
            (u_star*C + u_star/0.4.*log(y/delta)).*(y > delta);

    %高出的高度
    delta_h = Dp*sin(theta);

    %颗粒中心的流速
    u_center = u_y(0.25*Dp + delta_h);

    %积分上下限
    y1 = 0.25*Dp;
    y2 = 0.25*Dp + delta_h;

    %某一高度迎流面积函数
    len_y = @(y) 2*sqrt(Rp^2 - (y - delta_h + 0.75*Dp - Rp).^2);
    S = integral(len_y,y1,y2);

    %迎流面积的另外表示方式
    gamma = acos(2*sin(theta) - 1);
    A = Rp^2*(pi - gamma + sin(gamma)*cos(gamma));

    %速度面积积分的计算
    uA_y = @(y) u_y(y).*len_y(y);
    u2A_y = @(y) u_y(y).^2.*len_y(y);
    uA = integral(uA_y,y1,y2);
    u2A = integral(u2A_y,y1,y2);

    %速度的面积平均值
    u_uniform = uA/A;
    u2_uniform = sqrt(u2A/A);
    
    %u_uniform = 5.52*u_star;

    %计算拖曳力系数
    Rep = u_uniform*Dp/nv;
    Cd1 = 0;
    Cd2 = 0;
    Cd = 0;
    if Rep < 5
        Cd1 = 25/Rep;
    elseif Rep >= 5 && Rep <= 50000
        Cd1 = 0.55 + 37/Rep^1.2 - 3.5/Rep^0.9;
    else
        Cd1 = 0.25;
    end
    
    Cd2 = (0.63 + 4.8/sqrt(Rep))^2;

    Cd = Cd2;
    CD = Cd;

    %计算拖曳力
%     chi = 3.7 - 0.65*exp(-0.5*(1.5 - log10(Rep))^2);
%     cof = n^(-chi + 1);

    %不再适用拖曳力系数的经验公式，而采用平均的系数来表征
    %dCdR = -(0.63 + 4.8/Rep^0.5)*4.8/Rep^1.5;
    %alp = Rep/Cd*(dCdR);
    %F_drag0 = 0.5*Cd*rhof*A*u_uniform^2*(1 + (1 + 2*alp)*sigma_u^2/u_uniform^2);
    F_drag0 = 0.5*Cd*rhof*A*u_uniform^2;
    F_lift = 0.5*CL*rhof*A*(u_uniform^2 + sigma_u^2);  %时均升力比原来要大
    %F_lift = 0.5*CL*rhof*A*u_uniform^2;

    %扩展为湍流的判据
%     BR = pi*1.65*g*Dp^3/(3*A*(Cd*tan(theta) + CL));
%     BL = pi*1.65*g*Dp^3/(3*pi*Rp^2*CL);

end