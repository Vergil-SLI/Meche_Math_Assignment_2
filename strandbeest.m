function strandbeest()
    %clearing bc our stuff keeps breaking:
    clear;

    %Set Parameters
    %initialize leg_params structure
    leg_params = struct();
    %number of vertices in linkage
    leg_params.num_vertices = 7;
    %number of links in linkage
    leg_params.num_linkages = 10;
    %length of crank shaft
    leg_params.crank_length = 15.0;
    %fixed position coords of vertex 0 aka center of rotation of the crank
    leg_params.vertex_pos0 = [0;0];
    %fixed position coords of vertex 2 aka center of rotation of the 1
    %stationary vertex of the strandbeest leg
    leg_params.vertex_pos2 = [-38.0;-7.8];
    
    %matrix relating links to vertices
    leg_params.link_to_vertex_list = ...
           [1, 3;... %the 2 vertices link 1 connect
            3, 4;... %the 2 vertices link 2 connect
            2, 3;... %the 2 vertices link 3 connect
            2, 4;... %the 2 vertices link 4 connect
            4, 5;... %the 2 vertices link 5 connect
            2, 6;... %the 2 vertices link 6 connect
            1, 6;... %the 2 vertices link 7 connect
            5, 6;... %the 2 vertices link 8 connect
            5, 7;... %the 2 vertices link 9 connect
            6, 7]; %the 2 vertices link 10 connect
            
    leg_params.link_lengths = ...
           [50.0,... %link 1 length
            55.8,... %link 2 length
            41.5,... %link 3 length
            40.1,... %link 4 length
            39.4,... %link 5 length
            39.3,... %link 6 length
            61.9,... %link 7 length
            36.7,... %link 8 length
            65.7,... %link 9 length
            49.0]; %link 10 length
    
    vertex_guess_coords = [...
    [ 0; 50];... %vertex 1 [x,y] guess
    [ -50; 0];... %vertex 2 [x,y] guess
    [ -50; 50];... %vertex 3 [x,y] guess
    [-100; 0];... %vertex 4 [x,y] guess
    [-100; -50];... %vertex 5 [x,y] guess
    [ -50; -50];... %vertex 6 [x,y] guess
    [ -50; -100]... %vertex 7 [x,y] guess
    ];
    
    % concatenate vertex_guess_coords into a vertical column
    [numRows, numCols] = size(vertex_guess_coords);
    guess = zeros(numRows*numCols, 1);
    tracker = 1;
    for i = 1:numRows
        for j = 1:numCols
            guess(tracker) = vertex_guess_coords(i, j);
            tracker = tracker + 1;
        end
    end
    
    thetanum = 400; % higher number = smoother animation but longer time
    crank_final_angle = 12*pi; % how many rotation the animation does
    theta = linspace(0,crank_final_angle,thetanum);

    video_example(leg_params,guess,theta,thetanum);
end


function length_errors = link_length_error(vertex_coords, leg_params)
    %vertex_coords: a column vector containing the (x,y) coordinates of vertices (x1; y1; x2; y2...)   
    %length_errors: a column vector describing the current distance error of the ith link

    link_to_vertex_list = leg_params.link_to_vertex_list;
    link_lengths = leg_params.link_lengths;
    length_errors = zeros(leg_params.num_linkages, 1);

    for i = 1:leg_params.num_linkages
        vertices = link_to_vertex_list(i,:); % the 2 vertex connected by link i
        A = [vertex_coords(vertices(1)*2 - 1), vertex_coords(vertices(1)*2)]; % [XA, YA]
        B = [vertex_coords(vertices(2)*2 - 1), vertex_coords(vertices(2)*2)]; % (XB, YB)
        fi = (A(1)-B(1))^2 + (A(2)-B(2))^2 - link_lengths(i)^2;

        length_errors(i) = fi; 
    end
end


function coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta)
    %vertex_coords: a column vector containing the (x,y) coordinates of vertices (x1; y1; x2; y2...) 
    %coord_errors: a column vector of the diff. between the coords of vertices 1 & 2 vs what they should be
    vertex_1x = vertex_coords(1);
    vertex_1y = vertex_coords(2);
    vertex_2x = vertex_coords(3);
    vertex_2y = vertex_coords(4);

    vertex_1x_bar = leg_params.vertex_pos0(1) + leg_params.crank_length * cos(theta);
    vertex_1y_bar = leg_params.vertex_pos0(2) + leg_params.crank_length * sin(theta);
    vertex_2x_bar = leg_params.vertex_pos0(1) + leg_params.vertex_pos2(1);
    vertex_2y_bar = leg_params.vertex_pos0(2) + leg_params.vertex_pos2(2);

    coord_errors = [vertex_1x - vertex_1x_bar; vertex_1y - vertex_1y_bar;vertex_2x - vertex_2x_bar;vertex_2y - vertex_2y_bar];
end


function error_vec = linkage_error_func(vertex_coords, leg_params, theta)
    %concatenate length_errors and coord_errors
    distance_errors = link_length_error(vertex_coords, leg_params);
    coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta);
    error_vec = [distance_errors;coord_errors];
end


function [vertex_coords_root, exit_flag] = compute_coords(vertex_coords_guess, leg_params, theta)
    % using multi-newton to finally solve for them coordinates
    solver_params = struct();
    strandbeest_error = @(v) linkage_error_func(v, leg_params, theta);

    [vertex_coords_root, exit_flag] = multi_newton(strandbeest_error,vertex_coords_guess,solver_params);
end

% Plots the linkage or vertex of leg depending on what you want ig
function leg_drawing = initialize_leg_drawing(leg_params)
    leg_drawing = struct();
    leg_drawing.linkages = cell(leg_params.num_linkages,1);
    
    for linkage_index = 1:leg_params.num_linkages
        leg_drawing.linkages{linkage_index} = line([0,0],[0,0],'color','b','linewidth',2);
    end
    
    leg_drawing.crank = line([0,0],[0,0],'color','b','linewidth',1.5);
    leg_drawing.vertices = cell(leg_params.num_vertices,1);
    
    for vertex_index = 1:leg_params.num_vertices
        leg_drawing.vertices{vertex_index} = line([0],[0],'marker',...
        'o','markerfacecolor','r','markeredgecolor','r','markersize',8);
    end
end


%Plots leg given coordinates yay
function update_leg_drawing(complete_vertex_coords, leg_drawing, leg_params,theta,vertex_coords_root)
    axis equal
    axis([-150 20 -120 50])

    %iterate through each link, and plot each link
    for linkage_index = 1:leg_params.num_linkages
        %put the x & y coordinates of both ends of the link in line_x and y
        vertices = leg_params.link_to_vertex_list(linkage_index, :);
        line_x = [complete_vertex_coords(vertices(1)*2 - 1), complete_vertex_coords(vertices(2)*2 - 1)];
        line_y = [complete_vertex_coords(vertices(1)*2), complete_vertex_coords(vertices(2)*2)];
        set(leg_drawing.linkages{linkage_index},'xdata',line_x,'ydata',line_y);
    end

    %iterate through each vertex, and plot each vertex
    for vertex_index = 1:leg_params.num_vertices
        %put the x & y coordinates of the vertex into dot_x and y
        dot_x = complete_vertex_coords(vertex_index*2 - 1);
        dot_y = complete_vertex_coords(vertex_index*2);
        set(leg_drawing.vertices{vertex_index},'xdata',dot_x,'ydata',dot_y);
    end

    %crank_x and crank_y should both be two element arrays
    %containing the x and y coordinates of the line segment describing the crank
    crank_x = [leg_params.vertex_pos0(1,:), complete_vertex_coords(1)]; 
    crank_y = [leg_params.vertex_pos0(2,:), complete_vertex_coords(2)];
    set(leg_drawing.crank,'xdata',crank_x,'ydata',crank_y);
   
end


%Making the video
function video_example(leg_params,vertex_guess_coords,theta,thetanum)
    leg_drawing = initialize_leg_drawing(leg_params);
    hold off;
   
    % variables needed for Finite Diff and Linear Algebra velocity calc
    x7 = [];
    y7 = [];
    dv_x7 = [];
    dv_y7 = [];
    
    hold on;
    path_plot = plot(0,0,'g');
    tangent_plot = plot(0,0,'r');
    %iterate through theta 
    for theta_iter = 1:length(theta)     
        % calculate vertex coords for current theta and plot it 
        [vertex_coords_root, ~] = compute_coords(vertex_guess_coords, leg_params, theta(theta_iter));
        vertex_guess_coords = vertex_coords_root;


        % Finite diff data collection
        x7 = [x7; vertex_coords_root(13)];
        y7 = [y7; vertex_coords_root(14)];
        
        % Linear Algebra data collection
        dVdtheta = compute_velocities(vertex_coords_root, leg_params, theta(theta_iter));
        dv_x7 = [dv_x7; dVdtheta(13)];
        dv_y7 = [dv_y7; dVdtheta(14)];

        update_leg_drawing(vertex_coords_root, leg_drawing, leg_params,theta,vertex_coords_root)
        
        set(path_plot,'xdata',x7,'ydata',y7)
        set(tangent_plot,'xdata',[vertex_coords_root(13), vertex_coords_root(13)+dVdtheta(13)/2],...
            'ydata',[vertex_coords_root(14),vertex_coords_root(14)+dVdtheta(14)/2]);
        drawnow;    
    end
    % Finite Differences Calc
    dx7 = [];
    dy7 = [];

    for i = 1:length(x7)-1
        dx_dtheta = (x7(i+1)-x7(i)) / (theta(i+1)-theta(i));
        dy_dtheta = (y7(i+1)-y7(i)) / (theta(i+1)-theta(i));

        dx7 = [dx7; dx_dtheta];
        dy7 = [dy7; dy_dtheta];
    end

    figure();
    p_theta = linspace(0,(2*pi),thetanum);
    %Plot dxtip/dtheta
    subplot(2,1,1);
    plot(p_theta(1:length(p_theta)-1),dx7)
    hold on;
    plot(p_theta(1:length(p_theta)),dv_x7,'--')
    title('dxtip/dtheta Comparison')
    legend('Finite Differences', 'Linear Algebra')
    xlabel('Theta')
    ylabel('dxtip')
    axis([0 (2*pi) -50 50])
    hold off;
 
    % Plot dytip/dtheta
    subplot(2,1,2);
    plot(p_theta(1:length(p_theta)-1),dy7)
    hold on;
    plot(p_theta(1:length(p_theta)),dv_y7,'--')
    title('dytip/dtheta Comparison')
    legend('Finite Differences', 'Linear Algebra')
    xlabel('Theta')
    ylabel('dytip')
    axis([0 (2*pi) -50 50])
    hold off;
end

%theta: the current angle of the crank
%dVdtheta: a column vector containing the theta derivates of each vertex coord
function dVdtheta = compute_velocities(vertex_coords, leg_params, theta)
        wrapper = @(v) link_length_error(v,leg_params);
        Jacob_error_func = approximate_jacobian(wrapper,vertex_coords);
        
        add = eye(4,14);
        M = [add; Jacob_error_func];
        B = zeros(14,1);
        B(1,:) = leg_params.crank_length * -sin(theta);
        B(2,:) = leg_params.crank_length * cos(theta);
        
        dVdtheta = M\B;
    end
    
    function J = approximate_jacobian(fun,x)
        J = [];
        h = 1e-6;
    
        for j = 1:length(x)
            basis_j = zeros(length(x), 1);
            basis_j(j) = 1;
            column = (fun(x + h*basis_j) - fun(x - h*basis_j)) / (2*h);
            J = [J, column];
        end
    end