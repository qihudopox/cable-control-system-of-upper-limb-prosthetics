function [v,a,angle] = shoulder_flexion_rotation_model(step_t,invtime)
var = 0.0163;
mean = 0.0008;
gain = 55000;

t = 0:step_t:invtime;
v = gain*exp(-(t-mean).^2/(2*var)).*t.^3;
a = gain*exp(-(t-mean).^2/(2*var)).*t.^3.*...
    (-2*(t-mean)./(2*var))+3*gain*t.^2.*exp(-(t-mean).^2./(2*var));
angle = zeros(1,length(t));

for n = 2:length(t)
angle(n) = angle(n-1)+v(n)*step_t;
end

figure
subplot(3,1,1);
plot(t,v,'r');
title('flexion angular velocity')
ylabel('degree/s')
subplot(3,1,2);
plot(t,a,'b');
title('flexion angular acceleration')
ylabel('degree/s^2')
subplot(3,1,3);
plot(t,angle,'g');
title('flexion angle')
xlabel('time:s')
ylabel('degree')


end