clc; close all; clear all; % Clears all windows, figures, and workspaces

t_init = 0; % Initial time
t_fin = 1; % Final time

for j = 1:1000
    
    % Generates a random pair (x,y) within a unit triangle
    x_init = rand; 
    y_init = rand;
    while x_init + y_init > 1
        x_init = rand;
        y_init = rand;
    end

    h = 0.001; % Determines step size for numerical integration

    spaces = (t_fin - t_init)/h;  % Number of steps

    % Initialize arrays for time, x, and y values
    t_vals = zeros(1, spaces+1); 
    x_vals = zeros(1, spaces+1);
    y_vals = zeros(1, spaces+1);

    % Set initial values
    t_vals(1) = t_init; 
    x_vals(1) = x_init;
    y_vals(1) = y_init;

    % Define the differential equations
    x_dot = @(x,y,t) x.*(3.*x+4.*y-3);
    y_dot = @(x,y,t) y.*(3.*x+4.*y-4);

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
    xlabel('$x$','Interpreter','Latex');
    ylabel('$y$','Interpreter','Latex');
    a = [0 0 1];
    b = [1 0 0];
    labels = {'R', 'O', 'M'};
    text(a,b,labels,'VerticalAlignment','bottom','HorizontalAlignment','right', 'FontSize', 30,'Color','red')
    hold on
end

    