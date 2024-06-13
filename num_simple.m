clear;
clc;

%基本条件
g = 9.8;
rhos = 2650;
rhof = 1000;
Dp = 0.0025;
Rp = Dp/2;
cross_section = pi*Rp^2;
Vp = 4/3*pi*Rp^3;
Mp = Vp*rhos;
Gp = Vp*(rhos - rhof)*g;
nv = 1e-6;

%切应力
shi = 0.20;
tau = shi*1650*g*Dp;
%tau = 3.3205; %Pa
F0 = tau*cross_section;
u_star = sqrt(tau/rhof);
sigma_u = 2*u_star;
shields = tau/(rhos - rhof)/g/Dp;
Rstar = u_star*Dp*2/nv;
St = sqrt(tau)*Dp*sqrt(2650)*1000;
D_star = Dp*(1.65*g/nv^2)^(1/3);

%表征平均的拖曳力系数
C = 0.3907*(D_star - 14.5)^1.62*(St/10).^(-1.62);

%孔隙率
n = 0.4;

theta0list1 = 0.3:0.3:87.6;
theta0list2 = 87.7:0.1:89.9;
theta0list_deg = [theta0list1,theta0list2];
theta0list_rad = theta0list_deg*pi/180;
slopelist = [];
slope = 0*pi/180;

%计算拖曳力和升力
drag_list = [];
lift_list = [];
u0_list = [];
cof_list = [];
CD_list = [];

%初始角度不同，拖曳力也不同
for index = 1:length(theta0list_rad)
    [drag,lift,u0,CD] = calc_dragforce(Rp,tau,theta0list_rad(index),n);
    drag_list = [drag_list,drag];
    lift_list = [lift_list,lift];  %算时均升力
    u0_list = [u0_list,u0];  %计算得到的时均速度约为u*的5到8倍，随θ和D变化不大 5.52u*合理
    CD_list = [CD_list,CD];
end
u0max = max(u0_list);
u_ustar = u0_list/u_star;
Rep = u0_list*Dp/nv;

%给出theta0的分布
delta_deg1 = 0.3;
delta_deg2 = 0.1;


%根据组构
ac = 0.9;
thetac_deg = 0;
thetac = thetac_deg*pi/180;
pdfcn1 = 2*(1 + ac.*cos(2*(theta0list1*pi/180 - thetac)))/pi;
pdfcn2 = 2*(1 + ac.*cos(2*(theta0list2*pi/180 - thetac)))/pi;
%指数分布
lamda = 3;
pdf_exp1 = lamda.*exp(-lamda*theta0list1*pi/180);  %认为对称
pdf_exp2 = lamda.*exp(-lamda*theta0list2*pi/180);
scope_exp = sum(pdf_exp1.*delta_deg1)*pi/180 + sum(pdf_exp2.*delta_deg2)*pi/180;
mult = 1/scope_exp;
pdf_exp1 = pdf_exp1*mult;
pdf_exp2 = pdf_exp2*mult;
%假定暴露高度均匀分布 可推得角度分布为cosy
pdf_uni1 = cos(theta0list1*pi/180);
pdf_uni2 = cos(theta0list2*pi/180);
scope_uni = sum(pdf_uni1)*delta_deg1*pi/180 + sum(pdf_uni2)*delta_deg2*pi/180;
%对暴露高度施加一个线性分布 可推得角度分布
beta = 0.1;
pdf_buni1 = (1 + beta)*cos(theta0list1*pi/180) - beta*sin(2*theta0list1*pi/180);
pdf_buni2 = (1 + beta)*cos(theta0list2*pi/180) - beta*sin(2*theta0list2*pi/180);
scope_buni = sum(pdf_buni1)*delta_deg1*pi/180 + sum(pdf_buni2)*delta_deg2*pi/180;

% hold on
% plot(theta0list_deg,[pdfcn1,pdfcn2]);
% plot(theta0list_deg,[pdf_exp1,pdf_exp2]);
% plot(theta0list_deg,[pdf_uni1,pdf_uni2]);
% plot(theta0list_deg,[pdf_1uni1,pdf_1uni2]);

%定义要求的变量
EVi_list = [];  %没乘概率

%统计角度和速度
thetatot_list = [];
vy_list = [];
v_list = [];

L1 = length(theta0list1);
L2 = length(theta0list_deg);
flag = zeros(L2,1)';

%在主程序中进行循环计算侵蚀速率,只有theta0不同
%上浮1 滚动2 未拖动0 拖动但速度为负3
for i = 1:L2
    uu = u0_list(i);
    Fd = drag_list(i)*C;
    FL = lift_list(i);  %是时均升力，多出了速度波动的影响
    if Fd <= 0
        thetatot_list = [thetatot_list,0];
        vy_list = [vy_list,0];
        v_list = [v_list,0];
        EVi_list = [EVi_list,0];
        disp('拖曳力为负');
        continue;
    end
    bili = (Gp - FL)/(Fd);
    fai = atan(bili);
    
    if FL >= Gp
        disp([num2str(theta0list_deg(i)),'°直接上浮']);
        vy = sqrt(2*2*Rp*(FL - Gp)/Mp);
        thetatot_list = [thetatot_list,0];
        vy_list = [vy_list,vy];
        v_list = [v_list,vy];
        EVi_list = [EVi_list,vy*(1 - n)];
        
        flag(i) = 1;
        continue
    end
    
    if fai >= theta0list_rad(i)
        %disp([num2str(theta0list_deg(i)),'°无法拖动']);
        thetatot_list = [thetatot_list,0];
        vy_list = [vy_list,0];
        v_list = [v_list,0];
        EVi_list = [EVi_list,0];
        
        flag(i) = 0;
        continue
    end
    %二维
    co1 = Fd*cos(theta0list_rad(i)) + (Gp - FL)*sin(theta0list_rad(i));
    co2 = 2*Rp*10/17/Mp;
    v = sqrt(co1*co2);
    thetatot = real(acos(10/17*cos(theta0list_rad(i) - fai)) + fai);
    %三维有明显改善,但整体减小
%     co1 = Fd*cos(theta0list_rad(i)) + (Gp - FL)*sin(theta0list_rad(i));
%     co2 = sqrt(3)*Rp*30/53/Mp;
%     v = sqrt(co1*co2);
%     thetatot = acos(30/53*cos(theta0list_rad(i) - fai)) + fai;
    vy = real(v*cos(thetatot));
    
    thetatot_list = [thetatot_list,thetatot];
    vy_list = [vy_list,vy];
    v_list = [v_list,v];
    EVi_list = [EVi_list,vy*(1 - n)];
    if vy <= 0
        flag(i) = 3;
    else
        flag(i) = 2;
    end
end

flag_roll = (flag == 2);
flag_sus = (flag == 1);
flag_roll1 = flag_roll(1:L1);
flag_roll2 = flag_roll(1 + L1:end);
flag_sus1 = flag_sus(1:L1);
flag_sus2 = flag_sus(1 + L1:end);


EV_list = zeros(1,length(EVi_list));  %乘概率

% EV_list(1:L1) = EVi_list(1:L1).*pdfcn1*delta_deg1*pi/180;
% EV_list((L1+1):L2) = EVi_list((L1+1):L2).*pdfcn2*delta_deg2*pi/180;
% EV_list(1:L1) = EVi_list(1:L1).*pdf_exp1*delta_deg1*pi/180;
% EV_list((L1+1):L2) = EVi_list((L1+1):L2).*pdf_exp2*delta_deg2*pi/180;
EV_list(1:L1) = EVi_list(1:L1).*pdf_uni1*delta_deg1*pi/180;
EV_list((L1+1):L2) = EVi_list((L1+1):L2).*pdf_uni2*delta_deg2*pi/180;

EV_eft = find(EV_list > 0);

EV = 0;
if ~isempty(EV_eft)
    for j = 1:length(EV_eft)
        EV = EV + EV_list(EV_eft(j));  %单位面积单位时间侵蚀的体积
    end
else
    disp('竖向速度均为负');
    EV = 0;
end
Em = EV*rhos;
thetatot_list = thetatot_list*180/pi;

%输出结果
disp(Em);

% hold on;
% plot(theta0list_deg,vy_list);
% xlabel('旋转角（°）');
% ylabel('颗粒竖向速度（m/s）');