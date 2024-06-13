function [Em,Rstar,thetatot_list] = num_C(tau_sy,C_est,Dp)
    %基本条件
    g = 9.8;
    rhos = 2650;
    rhof = 1000;
    Rp = Dp/2;
    cross_section = pi*Rp^2;
    Vp = 4/3*pi*Rp^3;
    Mp = Vp*rhos;
    Gp = Vp*(rhos - rhof)*g;
    nv = 1e-6;

    %切应力
    tau = tau_sy;  %Pa
    F0 = tau*cross_section;
    u_star = sqrt(tau/rhof);
    sigma_u = 2*u_star;
    shields = tau/(rhos - rhof)/g/Rp/2;
    Rstar = u_star*2*Rp*2/nv;

    %表征平均的拖曳力系数
    C = C_est;
    gam = C_est;
    ud = 4/Rp*(sqrt(nv^2 + 0.0139*8*Rp^3*1.65*g) - nv);

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
        u0_list = [u0_list,u0]; %计算得到的时均速度约为u*的5到8倍，随θ和D变化不大 5.52u*合理
        CD_list = [CD_list,CD];
    end
    u0max = max(u0_list);

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
    lamda = 2;
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
            disp('直接上浮');
            vy = 0;
            thetatot_list = [thetatot_list,0];
            vy_list = [vy_list,vy];
            v_list = [v_list,vy];
            EVi_list = [EVi_list,vy*(1 - n)];
            
            flag(i) = 1;
            continue
        end

        if fai >= theta0list_rad(i)
            %disp('无法拖动');
            thetatot_list = [thetatot_list,0];
            vy_list = [vy_list,0];
            v_list = [v_list,0];
            EVi_list = [EVi_list,0];
            
            flag(i) = 0;
            continue
        end
        
        co1 = Fd*cos(theta0list_rad(i)) + (Gp - FL)*sin(theta0list_rad(i));
        co2 = 2*Rp*10/17/Mp;
        v = sqrt(co1*co2);
        thetatot = real(acos(10/17*cos(theta0list_rad(i) - fai)) + fai);
%         co1 = Fd*cos(theta0list_rad(i)) + (Gp - FL)*sin(theta0list_rad(i));
%         co2 = sqrt(3)*Rp*30/53/Mp;
%         v = sqrt(co1*co2);
%         thetatot = acos(30/53*cos(theta0list_rad(i) - fai)) + fai;
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

    
    %不同初始运动的分界点
    critical_roll = 0;
    critical_sus = 315;
    for k = 1:length(flag)
%         if flag(k) == 0
%             continue
%         elseif flag(k) == 3 || flag(k) == 2 || flag(k) == 1
%             critical_roll = k;
%             break
%         end
%         if flag(k) == 1
%             critical_sus = k;
%             break
%         else
%             continue
%         end
     end

%     roll_sus = theta0list_deg(critical_roll);
%     roll_sus = theta0list_deg(critical_sus);

    EV_list = zeros(1,length(EVi_list));  %乘概率

    % EV_list(1:L1) = EVi_list(1:L1).*pdfcn1*delta_deg1*pi/180;
    % EV_list((L1+1):L2) = EVi_list((L1+1):L2).*pdfcn2*delta_deg2*pi/180;
%     EV_list(1:L1) = EVi_list(1:L1).*pdf_exp1*delta_deg1*pi/180;
%     EV_list((L1+1):L2) = EVi_list((L1+1):L2).*pdf_exp2*delta_deg2*pi/180;
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
    %输出结果
    disp(Em);
    
    %计算不同运动形式的百分比 
%     static_s = 0;
%     roll_s = 0;
%     sus_s = 0;
%     for i = 1:L1
%         if flag(i) == 1
%             sus_s = sus_s + pdf_uni1(i)*delta_deg1*pi/180;
%         elseif flag(i) == 0
%             static_s = static_s + pdf_uni1(i)*delta_deg1*pi/180;
%         else
%             roll_s = roll_s + pdf_uni1(i)*delta_deg1*pi/180;
%         end
%     end
%     for i = (L1+1):L2
%         if flag(i) == 1
%             sus_s = sus_s + pdf_uni2(i - L1)*delta_deg2*pi/180;
%         elseif flag(i) == 0
%             static_s = static_s + pdf_uni2(i - L1)*delta_deg2*pi/180;
%         else
%             roll_s = roll_s + pdf_uni2(i - L1)*delta_deg2*pi/180;
%         end
%     end
    
end