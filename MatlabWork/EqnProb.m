clear; clc;

R = 300; F = 0.466*R;

syms p ZA ZB;
eqn1 = sqrt(150^2 + (ZB + R - F)^2) == ZB - 2*ZA + F - R;
eqn2 = (150*300-2*p*ZB) / sqrt(150^2+ZB^2) == ...
        (150*300-2*p*(ZB+R-F)) / sqrt(150^2+(ZB+R-F)^2);
eqn3 = 2*p*ZB - 2*p*ZA == 150^2;

S = solve([eqn1 eqn2 eqn3], [p ZA ZB]);

disp("p:");
disp(double(S.p));

disp("ZA:");
disp(double(S.ZA));

disp("ZB:");
disp(double(S.ZB));

% disp([double(S.p), double(S.ZA), double(S.ZB)]);