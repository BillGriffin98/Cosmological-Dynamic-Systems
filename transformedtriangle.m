clc; close all; clear all; % Clears all windows, figures, and workspaces

t_init = 0; % Initial time
t_fin = 1; % Final time

for j = 1:100
    
    % Generates (x,y) such that x is in (0, 2) and y is in (0, 1)
    x_init = 2*rand;
    y_init = rand;
    while ((1/2).*((-2)+(-1).*x_init+y_init+(8.*y_init+((-2)+(-1).*x_init+y_init).^2).^(1/2)) + (1/2).*(x_init+(-1).*y_init+(1/2).*((-2)+(-1).*x_init+y_init+(8.*y_init+((-2)+(-1).*x_init+y_init).^2).^(1/2)))) > 1 
        x_init = rand;
        y_init = rand;
    end

    h = 0.001; % Determines step size for numerical integration

    spaces = (t_fin - t_init)/h; % Number of steps

    % Initialize arrays for time, x, and y values
    t_vals = zeros(1, spaces+1);
    x_vals = zeros(1, spaces+1);
    y_vals = zeros(1, spaces+1);

    % Set initial values
    t_vals(1) = t_init;
    x_vals(1) = x_init;
    y_vals(1) = y_init;

    % Define the differential equations
    x_dot = @(m,r,t) (1/2).*(2+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))).*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))).*((-4)+(3/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))+2.*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))))+(1/4).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2)).*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))).*((-3)+(3/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))+2.*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))));
    y_dot = @(m,r,t) (1/4).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2)).*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))).*((-4)+(3/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))+2.*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))))+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2)).*(1+(1/2).*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2)))).*((-3)+(3/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))+2.*(m+(-1).*r+(1/2).*((-2)+(-1).*m+r+(8.*r+((-2)+(-1).*m+r).^2).^(1/2))));

    % Runge-Kutta method for numerical integration
    for i = 1:spaces
            t_vals(i+1) = t_vals(i) + h;

            k1X = h*x_dot(x_vals(i), y_vals(i), t_vals(i));
            k1Y = h*y_dot(x_vals(i), y_vals(i), t_vals(i));

            k2X = h*x_dot(x_vals(i) + k1X/2, y_vals(i) + k1Y/2, t_vals(i) + h/2);
            k2Y = h*y_dot(x_vals(i) + k1X/2, y_vals(i) + k1Y/2, t_vals(i) + h/2);

            k3X = h*x_dot(x_vals(i) + k2X/2, y_vals(i) + k2Y/2, t_vals(i) + h/2);
            k3Y = h*y_dot(x_vals(i) + k2X/2, y_vals(i) + k2Y/2, t_vals(i) + h/2);

            k4X = h*x_dot(x_vals(i) + k3X, y_vals(i) + k3Y, t_vals(i) + h);
            k4Y = h*y_dot(x_vals(i) + k3X, y_vals(i) + k3Y, t_vals(i) + h);

            x_vals(i+1) = x_vals(i) + (1/6)*(k1X + 2*k2X + 2*k3X + k4X);
            y_vals(i+1) = y_vals(i) + (1/6)*(k1Y + 2*k2Y + 2*k3Y + k4Y);
        end

    % Plot the trajectory in the x-y plane for each iteration
    figure(1)
    plot(x_vals, y_vals);

    % Label the axis
    axis equal;
    set(gca,'FontSize',20)
    xlabel('$m$','Interpreter','Latex');
    ylabel('$r$','Interpreter','Latex');
    hold on
end

    