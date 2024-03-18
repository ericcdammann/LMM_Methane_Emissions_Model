function [value, isterminal, direction] = behavior(~, y)

value      = [y(3) >= 1.4 || y(3)>=y(2); y(1)>=0.99*y(4); y(1)<=0.01*y(4) && y(2) <= 1.4]; % Drown; erode; fill
isterminal = [1; 1; 1];
direction  = [1; 1; 1];