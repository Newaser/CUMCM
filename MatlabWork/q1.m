clc; clear;
syms x 
R = 300;
F = 0.466*R;
h = -0.6:0.01:0.6;
% f = h + R;
% k = f ./ h;
k = h ./(R - F);
S = int((abs((x^2-R-h)./(4*R.*k)) + sqrt(R^2-x^2)), x, [-150 150]);
plot(h, S);
xlabel('h');
ylabel('S');

[S0, idx_h0] = min(S);
h0 = -0.6 + 0.01*(idx_h0-1);
disp(double(S0));
fprintf("\n%f", S0);
disp(h0);