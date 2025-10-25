clc; clear all; close all;

w = 1;
t = 0:0.1:20; 
num = w^2;

% Input multiple zeta values
zeta_values = input('Enter zeta values as a vector (e.g. [0.1 0.5 1]): ');

figure; hold on;
for i = 1:length(zeta_values)
    z = zeta_values(i);
    den = [1 2*z*w w^2];
    
    sys = tf(num, den);           % Create transfer function
    
    y = step(sys, t);             % Get step response
    
    plot(t, y, 'DisplayName', ['\zeta = ' num2str(z)]);
    
    info = stepinfo(y, t);        % Use y and t to get step info
    fprintf('Zeta = %.2f, RiseTime = %.3f, SettlingTime = %.3f\n', z, info.RiseTime, info.SettlingTime);
end
hold off;

legend show;
xlabel('Time (s)');
ylabel('Response');
title('Step Response for Multiple Zeta Values');
grid on;
