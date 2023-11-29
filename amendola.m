clc; close all; clear all; % Clears all windows, figures, and workspaces

t_init = 0; % Initial time
t_fin = 1; % Final time

for j = 1:2000

    % Generates (x,y) such that x is in (-1.5, 5) and y is in (-5, 11)
    y_init = 3-4.5*rand;
    x_init = 3-8*rand;
   
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
    x_dot = @(m,r,t) -4 + 2.*r + 3.*m + m.^2 - m.*r;
    y_dot = @(m,r,t) -(m.*r)-2.*r.*(r-2);

    % Runge-Kutta method for numerical integration
    for i = 1:spaces
        t_vals(i+1) = t_vals(i) + h;

        k1X = h*x_dot(x_vals(i), y_vals(i), t_vals(i));
        k1Y = h*y_dot(x_vals(i), y_vals(i), t_vals(i));

        k2X = h*x_dot(x_vals(i) + k1X/2, y_vals(i) + k1X/2, t_vals(i) + h/2);
        k2Y = h*y_dot(x_vals(i) + k1X/2, y_vals(i) + k1X/2, t_vals(i) + h/2);

        k3X = h*x_dot(x_vals(i) + k2X/2, y_vals(i) + k2X/2, t_vals(i) + h/2);
        k3Y = h*y_dot(x_vals(i) + k2X/2, y_vals(i) + k2X/2, t_vals(i) + h/2);

        k4X = h*x_dot(x_vals(i) + k3X, y_vals(i) + k3X, t_vals(i) + h);
        k4Y = h*y_dot(x_vals(i) + k3X, y_vals(i) + k3X, t_vals(i) + h);


        x_vals(i+1) = x_vals(i) + (1/6)*(k1X + 2*k2X + 2*k3X + k4X);
        y_vals(i+1) = y_vals(i) + (1/6)*(k1Y + 2*k2Y + 2*k3Y + k4Y);
        
    end

    % Plot the trajectory in the y-x plane for each iteration
    figure(1) 
    plot(x_vals, y_vals); 
    hold on; 
   
    % Label axis
    axis( [ -5 3 -1.5 3 ] ); 
    set(gca,'FontSize',20)
    xlabel('$x_1$','Interpreter','Latex');
    ylabel('$x_3$','Interpreter','Latex');
    
        
    

end
