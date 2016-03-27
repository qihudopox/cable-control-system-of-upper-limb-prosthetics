%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on constrains, calculate the values of c,d,e,f
%
% Assuming cable is stiff and will not deform(lengthen) when stretch
%
% Analyzing equilibrium status(no transient analysis)
%
% Model the case of lifting the upper arm from about vertical
% to about level when holding things with mass wlo
% obtain the maximum force through the process
%
% Model the case of opening the finger
%
% Model the dynamic respond
% NB:shoulder flexion angular velocity and acceleration is coded in
% funcion shoulder_flexion_rotation_model. This model can be modified
% within the function
%
% All output units are in S.I. units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system choice
geometry = 1; %1:calculate geometry; 0: skip the step
equilibrium_fel = 1; %1:calculate equilibrium fel; 0: skip the step
equilibrium_fth = 1; %1:calculate equilibrium fth; 0: skip the step
dynamic = 1;  %1:activate dynamic model; 0: skip the step
friction = 1; %1: rotational friction depends on angular velocity
%2: rotational friction is a constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%determined independent variables
%see logbook page39, 41, 49, 51 for meanings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit:cm
ext = 7.5;
thl = 8;
ul = 18;
ll = 13;
a = 4;
b = 4;
p = 4; % p < thl
cg0 = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit = degree
thang = 10;
ulang1 = 30;
ulang2 = 45.5;
llangmax = 170;
llangmin = 130;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit:kg
mlh = 0.2;
mlo = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unit:N
firb = 1.7; % it is chosen such that force applied
% to open hand is greater than force applied to lift arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit:N/cm
k = 3;  % it is chosen such that force applied
% to open hand is greater than force applied to lift arm
kcbel = 150;    %cable elastic constant = E * cable_length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unit: kg*m^2
ilh = 0.0042;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unit: N*s
mfri = 0.3;   %angular friction coefficient, affected by rotational speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unit: N*m
elbow_frictional_moment = 0.3; %elbow friction induced anti torque
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step
step_for_angle = 1; % unit:degree
step_for_time = 0.001; %unit:second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%unit:s
invtime = 3; %amount of time for upper arm to achieve target position
dynamic_tracktime = 3;  %length of time to track the dynamic response of
%lower arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do not edit anything below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%independent determined modify unit
wlh = mlh*9.81;
wlo = mlo*9.81;
step_for_angle = step_for_angle/180*pi;
thang = thang/180*pi;
ulang1 = ulang1/180*pi;
ulang2 = ulang2/180*pi;
llangmax = llangmax/180*pi;
llangmin = llangmin/180*pi;
no_of_gamma = ulang1/step_for_angle+1;
ext = ext/100;
thl = thl/100;
ul = ul/100;
ll = ll/100;
a = a/100;
b = b/100;
p = p/100; % p < thl
cg0 = cg0/100;
k = k*100;
kcbel = kcbel*100;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solving theta from ext
if (geometry == 1)
    theta = acos(1-(ext/(sqrt(2)*thl))^2)+thang;
    
    for count1 = 1:length(theta)
        if (theta(count1)>=0)
            sstheta = theta(count1);
        end
    end
    
    theta = sstheta/pi*180;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %solving c,d,e,f from range of motion
    syms d
    c = ul-b;
    [sd] = solve(sqrt(c^2+d^2-2*cos(llangmax)*c*d)-...
        sqrt(c^2+d^2-2*cos(llangmin)*c*d)==...
        sqrt(b^2+a^2-2*cos(ulang1)*a*b)-abs(b-a),d);
    d=double(sd);
    e = ll-d;
    
    syms f
    
    [sf]=solve(sqrt(e^2+f^2-2*cos(pi-thang)*e*f)-...
        sqrt(e^2+f^2-2*cos(pi-sstheta)*e*f)==...
        sqrt(b^2+a^2-2*cos(ulang2)*a*b)-sqrt(b^2+a^2-2*cos(ulang1)*a*b));
    
    f = double(sf);
    
    display(a);
    display(b);
    display(c);
    display(d);
    display(e);
    display(f);
    display(theta);
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving max force and its respective angle for the motion
% of lifting lower limb with a load mass of wlo
if (equilibrium_fel == 1)
    
    for count2 = 1:no_of_gamma
        
        gamma(count2) = (count2-1)*step_for_angle;
        
        
        ssalpha = acos(-((abs(b-a)-sqrt(b^2+a^2-2*a*b*cos(gamma(count2)))+...
            sqrt(c^2+d^2-2*c*d*cos(llangmax)))^2-(c^2+d^2))/(2*c*d));
        
        
        for count3 = 1:length(ssalpha)
            if ssalpha(count3)>=llangmin-0.1/180*pi
                if ssalpha(count3)<= llangmax+0.1/180*pi
                    sssalpha = ssalpha(count3);
                end
            end
        end
        
        solalpha(count2) = sssalpha;
        
        solbeta(count2) = pi-solalpha(count2)-...
            asin(sin(solalpha(count2))*d/...
            sqrt(d^2+c^2-2*d*c*cos(solalpha(count2))));
        
        
        fel = (wlh*cos(solalpha(count2)-pi/2-gamma(count2))*cg0+...
            wlo*cos(solalpha(count2)-pi/2-gamma(count2))*(thl+ll))/...
            (sin(solbeta(count2))*d);
        
        
        solfel(count2) = fel;
        
        
        
    end
    
    
    maxfel = solfel(1);
    
    for count4 = 1:no_of_gamma
        if solfel(count4) >= maxfel
            maxfel = solfel(count4);
            cn1 = count4;
        end
    end
    
    max_force_with_gamma = ((cn1-1)*step_for_angle)/pi*180;
    max_force_with_alpha = solalpha(cn1)/pi*180;
    max_force_with_beta = solbeta(cn1)/pi*180;
    
    figure
    plot(gamma/pi*180,solfel,'b');
    title('static cable tensile force to lift forearm');
    xlabel('upper arm rotation : degree');
    ylabel('static cable tensile force : N');
    hold on
    plot(max_force_with_gamma,maxfel,'r+');
    
    
    display(max_force_with_gamma);
    display(max_force_with_alpha);
    display(max_force_with_beta);
    display(maxfel);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solving max force to open thumb

if (equilibrium_fth == 1)
    maxthetacount = ceil((sstheta-thang)/step_for_angle);
    for count5 = 1: maxthetacount
        mtheta(count5) = ((count5-1)*step_for_angle+thang);
        mpsi = asin(f*sin(pi-mtheta(count5))/...
            sqrt(f^2+e^2-2*e*f*cos(pi-mtheta(count5))));
        extrb = sqrt(p^2+p^2-2*p*p*cos(mtheta(count5)))-...
            sqrt(p^2+p^2-2*p*p*cos(thang));
        phi = (pi-mtheta(count5))/2;
        fth(count5) = (firb+k*extrb)*f*sin(phi)/...
            (sin(mpsi)*e);
        frb(count5) = firb+k*extrb;
    end
    maxfth = fth(1);
    for count6 = 1:maxthetacount
        if fth(count6) >= maxfth
            maxfth = fth(count6);
            cn2 = count6;
        end
    end
    
    max_force_with_mtheta = ((cn2-1)*step_for_angle+thang)/pi*180;
    
    figure
    plot(mtheta/pi*180,fth,'b');
    hold on
    plot(mtheta/pi*180,frb,'r')
    title('static cable tensile force to open hand');
    xlabel('thumb rotation : degree');
    ylabel('force : N');
    hold on
    plot(max_force_with_mtheta,maxfth,'r+');
    
    display(maxfth);
    display(max_force_with_mtheta);
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model for moving system
if dynamic == 1
    
    %calculate cg including load
    cg = (cg0 *mlh + (ll+thl)*mlo)/(mlh+mlo);
    
    % numerical modeling of input shoulder flexion rotation
    % angular velocity and angular acceleration
    [omegaulm,alphaulm,thetaulm] =...
        shoulder_flexion_rotation_model(step_for_time,invtime);
    
    omegaf = 0;
    thetaf = llangmax;
    
    syms alphaf
    
    for count7 = 1:length(omegaulm)
        omegaul = omegaulm(count7)/180*pi;
        thetaul = thetaulm(count7)/180*pi;
        alphaul = alphaulm(count7)/180*pi;
        
        thetaf(count7+1) = thetaf(count7) +...
            (omegaul-omegaf(count7))*step_for_time;
        
        
        extcbel(count7) = sqrt(a^2+b^2-2*a*b*cos(thetaul))+...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count7+1)))-abs(b-a)-...
            sqrt(c^2+d^2-2*c*d*cos(llangmax));
        
        sirat = sin(thetaf(count7))/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count7)));
        
        if friction ==1
            fri = mfri*(omegaf(count7)-omegaul);
        elseif friction ==2
            if abs(omegaf(count7)-omegaul) == 0
                fri = 0;
            else
                fri = elbow_frictional_moment*(omegaf(count7)-omegaul)/...
                    abs((omegaf(count7)-omegaul));
            end
        end
        
        
        fulel =(kcbel*extcbel(count7)*c*(d-cg)*sin(thetaf(count7+1))/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count7+1)))-fri-...
            (ilh +mlh*(cg-cg0)^2+mlo*(ll+thl-cg)^2)*alphaf)/...
            (sin(pi-thetaf(count7+1))*cg);
        
        
        [salphaf] = solve((mlh+mlo)*(omegaul^2*(b+c)*d*sirat+...
            alphaul*(b+c)*cos(asin(d*sirat))-...
            omegaf(count7)^2*cg*(b+c)*sirat+...
            alphaf*cg*cos(asin((b+c)*sirat)))==fulel*d*sirat-...
            kcbel*extcbel(count7)*b*d*sirat/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count7+1)))-...
            (wlh+wlo)*sin(thetaul+asin(d*sirat)),alphaf);
        salphaf = double(salphaf);
        
        if thetaf(count7+1)>llangmax
            thetaf(count7+1) = llangmax;
            ssalphaf(count7)= alphaul;
            omegaf(count7+1) = omegaul;
        else
            ssalphaf(count7) = salphaf;
            omegaf(count7+1) = omegaf(count7) + ssalphaf(count7)*...
                step_for_time;
        end
        
        %         ssfulel(count7) = ((mlh+mlo)*(omegaul^2*(b+c)*d*sirat+...
        %             alphaul*(b+c)*cos(asin(d*sirat))-...
        %             omegaf(count7)^2*cg*(b+c)*sirat+...
        %             ssalphaf(count7)*cg*cos(asin((b+c)*sirat))) + ...
        %             kcbel*extcbel(count7)*b*d*sirat/...
        %             sqrt(c^2+d^2-2*c*d*cos(thetaf(count7+1)))+...
        %             (wlh+wlo)*sin(thetaul+asin(d*sirat)))/d*sirat;
        
        ssfulel(count7) =(kcbel*extcbel(count7)*c*(d-cg)*sin(thetaf(count7+1))/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count7+1)))-fri-...
            (ilh +mlh*(cg-cg0)^2+mlo*(ll+thl-cg)^2)*ssalphaf(count7))/...
            (sin(pi-thetaf(count7+1))*cg);
        
        
        
    end
    
    
    
    for count8 = ...
            (invtime/step_for_time+1):(dynamic_tracktime/step_for_time+1)
        
        
        thetaf(count8+1) = thetaf(count8) +...
            (omegaul-omegaf(count8))*step_for_time;
        
        
        extcbel(count8) = sqrt(a^2+b^2-2*a*b*cos(thetaul))+...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count8+1)))-abs(b-a)-...
            sqrt(c^2+d^2-2*c*d*cos(llangmax));
        
        sirat = sin(thetaf(count8))/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count8)));
        
        if friction ==1
            fri = mfri*(omegaf(count8)-omegaul);
        elseif friction ==2
            if abs(omegaf(count8)-omegaul) ==0
                fri = 0;
            else
                fri = elbow_frictional_moment*(omegaf(count8)-omegaul)/...
                    abs((omegaf(count8)-omegaul));
            end
        end
        
        
        fulel =(kcbel*extcbel(count8)*c*(d-cg)*sin(thetaf(count8+1))/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count8+1)))-fri-...
            (ilh +mlh*(cg-cg0)^2+mlo*(ll+thl-cg)^2)*alphaf)/...
            (sin(pi-thetaf(count8+1))*cg);
        
        
        [salphaf] = solve((mlh+mlo)*(omegaul^2*(b+c)*d*sirat+...
            alphaul*(b+c)*cos(asin(d*sirat))-...
            omegaf(count8)^2*cg*(b+c)*sirat+...
            alphaf*cg*cos(asin((b+c)*sirat)))==fulel*d*sirat-...
            kcbel*extcbel(count8)*b*d*sirat/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count8+1)))-...
            (wlh+wlo)*sin(thetaul+asin(d*sirat)),alphaf);
        salphaf = double(salphaf);
        
        if thetaf(count8+1)>llangmax
            thetaf(count8+1) = llangmax;
            ssalphaf(count8)= alphaul;
            omegaf(count8+1) = omegaul;
        else
            ssalphaf(count8) = salphaf;
            omegaf(count8+1) = omegaf(count8) + ssalphaf(count8)*...
                step_for_time;
        end
        
        %         ssfulel(count8) = ((mlh+mlo)*(omegaul^2*(b+c)*d*sirat+...
        %             alphaul*(b+c)*cos(asin(d*sirat))-...
        %             omegaf(count8)^2*cg*(b+c)*sirat+...
        %             ssalphaf(count8)*cg*cos(asin((b+c)*sirat))) + ...
        %             kcbel*extcbel(count8)*b*d*sirat/...
        %             sqrt(c^2+d^2-2*c*d*cos(thetaf(count8+1)))+...
        %             (wlh+wlo)*sin(thetaul+asin(d*sirat)))/d*sirat;
        
        ssfulel(count8) =(kcbel*extcbel(count8)*c*(d-cg)*sin(thetaf(count8+1))/...
            sqrt(c^2+d^2-2*c*d*cos(thetaf(count8+1)))-fri-...
            (ilh +mlh*(cg-cg0)^2+mlo*(ll+thl-cg)^2)*ssalphaf(count8))/...
            (sin(pi-thetaf(count8+1))*cg);
        
        omegaulm(count8+1)=0;
        thetaulm(count8+1)=thetaulm(count8);
        alphaulm(count8+1)=0;
        
    end
    
    
    
    t = 0:step_for_time:dynamic_tracktime;
    omegaf = omegaf.*180./pi;
    ssalphaf = ssalphaf.*180./pi;
    thetaf = thetaf.*180./pi;
    
    figure
    subplot(3,1,1);
    plot(t,omegaf(1:(length(omegaf)-1)),'r');
    hold on
    plot(t,omegaulm(1:(length(omegaf)-1)),'c');
    if friction ==1
        titlstr = sprintf('cable constant = %g N/m   rotational friction coefficient = %g N*s d = %g cm \n flexion angular velocity',kcbel,mfri,d);
    elseif friction ==2
        titlstr = sprintf('cable constant = %g N/m   rotational friction = %g N*m d = %g cm \n flexion angular velocity',kcbel,elbow_frictional_moment,d);
    end
    title(titlstr);
    ylabel('degree/s');
    subplot(3,1,2);
    plot(t,ssalphaf,'b');
    hold on
    plot(t,alphaulm(1:(length(omegaf)-1)),'c');
    ylabel('degree/s^2');
    title('flexion angular acceleration');
    subplot(3,1,3);
    plot(t,thetaf(2:length(thetaf)),'g');
    hold on
    plot(t,thetaulm(1:(length(omegaf)-1)),'c');
    title('flexion angle');
    ylabel('degree');
    xlabel('time: s');
    
    figure
    plot(t,extcbel.*kcbel);
    title('cable tensile force');
    ylabel('N');
    xlabel('time: s');
    
    figure
    plot(t,ssfulel);
    title('elbow force');
    ylabel('N');
    xlabel('time: s')
    
    
end

