function projectile()
    solver_params = struct();
    solver_params.numerical_diff = 1;
    
    [x_guess, exit_flag] = multi_newton(@object_dist, [pi/4; 10], solver_params)
    
    run_simulation(x_guess(1),x_guess(2)) % theta, tc
end

function fval = object_dist(x) 
    % x = theta, t
    fval = projectile_traj(x(1),x(2)) - target_traj(x(2));
end

%projectile motion function
%theta is angle projectile is fired at (in radians)
%t is time in seconds
function V_p = projectile_traj(theta,t)
    g = 2.3; %gravity in m/sˆ2
    v0 = 14; %initial speed in m/s
    px0 = 2; %initial x position
    py0 = 4; %initial y position
    %compute position vector
    V_p = [v0*cos(theta)*t+px0; -.5*g*t.^2+v0*sin(theta)*t+py0];
end

%function describing motion of a flying target
%t is time in seconds
function V_t = target_traj(t)
    a1 = 7; %x amplitude in meters
    b1 = 9; %y amplitude meters
    omega1 = 3; %frequency in rad/sec
    phi1 = -pi/7; %phase shift in radians
    a2 = 2; %x amplitude in meters
    b2 = .7; %y amplitude meters
    omega2 = 5; %frequency in rad/sec
    phi2 = 1.5*pi; %phase shift in radians
    x0 = 28; %x offset in meters
    y0 = 21; %y offset in meters
    %compute position vector
    V_t = [a1*cos(omega1*t+phi1)+a2*cos(omega2*t+phi2)+x0;...
    b1*sin(omega1*t+phi1)+b2*sin(omega2*t+phi2)+y0];
end

%visualize the motion of the projectile and the target
%theta is the initial firing angle of the projectile
%t_c is the predicted collision time
function run_simulation(theta,t_c)
    %create the plot window, set the axes size, and add labels
    hold on;
    axis equal; axis square;
    axis([0,50,0,50])
    xlabel('x (m)')
    ylabel('y (m)')
    title('Simulation of Projectile Shot at Target')
    %initialize plots for the projectile/target and their paths
    traj_line_proj = plot(0,0,'g--','linewidth',2);
    traj_line_targ = plot(0,0,'k--','linewidth',2);
    t_plot = plot(0,0,'bo','markerfacecolor','b','markersize',8);
    p_plot = plot(0,0,'ro','markerfacecolor','r','markersize',8);
    %position lists
    %used for plotting paths of projectile/target
    V_list_proj = [];
    V_list_targ = [];
    %iterate through time until a little after the collision occurs
    for t = 0:.005:t_c+1.5
        %set time so that things freeze once collision occurs
        t_input = min(t,t_c);
        %compute position of projectile and target
        V_p = projectile_traj(theta,t_input);
        V_t = target_traj(t_input);
        %update the position lists
        V_list_proj(:,end+1) = V_p;
        V_list_targ(:,end+1) = V_t;
        %index used for tail of target path
        i = max(1,size(V_list_targ,2)-300);
        %update plots
        set(t_plot,'xdata',V_t(1),'ydata',V_t(2));
        set(p_plot,'xdata',V_p(1),'ydata',V_p(2));
        set(traj_line_proj,'xdata',V_list_proj(1,:),'ydata',V_list_proj(2,:));
        set(traj_line_targ,'xdata',V_list_targ(1,i:end),'ydata',V_list_targ(2,i:end));
        %show updated plots
        drawnow;
    end
end