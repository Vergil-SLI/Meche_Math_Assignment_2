function strandbeest()
    %Set Parameters
    %initialize leg_params structure
    leg_params = struct();
    %number of vertices in linkage
    leg_params.num_vertices = 7;
    %number of links in linkage
    leg_params.num_linkages = 10;
    %length of crank shaft
    leg_params.crank_length = 15.0;
    %fixed position coords of vertex 0
    leg_params.vertex_pos0 = [0;0];
    %fixed position coords of vertex 2 ***********************
    leg_params.vertex_pos2 = [-38.0;-7.8];
    
    %matrix relating links to vertices
    leg_params.link_to_vertex_list = ...
           [1, 3;... %link 1 adjacency
            3, 4;... %link 2 adjacency
            2, 3;... %link 3 adjacency
            2, 4;... %link 4 adjacency
            4, 5;... %link 5 adjacency
            2, 6;... %link 6 adjacency
            1, 6;... %link 7 adjacency
            5, 6;... %link 8 adjacency
            5, 7;... %link 9 adjacency
            6, 7]; %link 10 adjacency
            
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
    [ 0; 50];... %vertex 1 guess
    [ -50; 0];... %vertex 2 guess
    [ -50; 50];... %vertex 3 guess
    [-100; 0];... %vertex 4 guess
    [-100; -50];... %vertex 5 guess
    [ -50; -50];... %vertex 6 guess
    [ -50; -100]... %vertex 7 guess
    ];
    
    theta = 0;
    
    % linkage_error_func(guess, leg_params, theta)
    [vertex_coords_root, ~] = compute_coords(vertex_guess_coords, leg_params, theta)

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

    coord_errors = [
        vertex_1x - vertex_1x_bar;...
        vertex_1y - vertex_1y_bar;...
        vertex_2x - vertex_2x_bar;...
        vertex_2y - vertex_2y_bar];
end


function error_vec = linkage_error_func(vertex_coords, leg_params, theta)
    %concatenate length_errors and coord_errors
    distance_errors = link_length_error(vertex_coords, leg_params);
    coord_errors = fixed_coord_error_func(vertex_coords, leg_params, theta);
    error_vec = [distance_errors;coord_errors];
end


function [vertex_coords_root, exit_flag] = compute_coords(vertex_coords_guess, leg_params, theta)
    % just find a config of vertices that fits the length bro
    solver_params = struct();
    strandbeest_error = @(v) linkage_error_func(v, leg_params, theta);
    
    [numRows, numCols] = size(vertex_coords_guess);
    guess = zeros(numRows*numCols, 1);
    tracker = 1;
    for i = 1:numRows
        for j = 1:numCols
            guess(tracker) = vertex_coords_guess(i, j);
            tracker = tracker + 1;
        end
    end

    [vertex_coords_root, exit_flag] = multi_newton(strandbeest_error,guess,solver_params);
end

%% EVERYTHING UNDERNEATH I HAVE NOT DONE, WILL DO TMR
%leg_drawing: a struct containing all the plotting objects for the linkage
% leg_drawing.linkages is a cell array, where each element corresponds
% to a plot of a single link (excluding the crank)
% leg_drawing.crank is a plot of the crank link
% leg_drawing.vertices is a cell array, where each element corresponds
% to a plot of one of the vertices in the linkage
function leg_drawing = initialize_leg_drawing(leg_params)
    leg_drawing = struct();
    leg_drawing.linkages = cell(leg_params.num_linkages,1);
    
    for linkage_index = 1:leg_params.num_linkages
        leg_drawing.linkages{linkage_index} = line([0,0],[0,0],'color','k','linewidth',2);
    end
    
    leg_drawing.crank = line([0,0],[0,0],'color','k','linewidth',1.5);
    leg_drawing.vertices = cell(leg_params.num_vertices,1);
    
    for vertex_index = 1:leg_params.num_vertices
        leg_drawing.vertices{vertex_index} = line([0],[0],'marker',...
        'o','markerfacecolor','r','markeredgecolor','r','markersize',8);
    end
end


%Updates the plot objects that visualize the leg linkage
%for the current leg configuration
%INPUTS:
%vertex_coords_root: a column vector containing the (x,y) coordinates of every vertex
%leg_drawing: a struct containing all the plotting objects for the linkage
% leg_drawing.linkages is a cell array, where each element corresponds
% to a plot of a single link (excluding the crank)
% leg_drawing.crank is a plot of the crank link
% leg_drawing.vertices is a cell array, where each element corresponds
% to a plot of one of the vertices in the linkage
% function update_leg_drawing(complete_vertex_coords, leg_drawing, leg_params)
%     %iterate through each link, and update corresponding link plot
%     for linkage_index = 1:leg_params.num_linkages
%         %linkage_index is the label of the current link
%         %your code here
%         %line_x and line_y should both be two element arrays containing
%         %the x and y coordinates of the line segment describing the current link
%         line_x = %your code here
%         line_y = %your code here
%         set(leg_drawing.linkages{linkage_index},'xdata',line_x,'ydata',line_y);
%     end
%     %iterate through each vertex, and update corresponding vertex plot
%     for vertex_index = 1:leg_params.num_vertices
%         %vertex_index is the label of the current vertex
%         %your code here
%         %dot_x and dot_y should both be scalars
%         %specifically the x and y coordinates of the corresponding vertex
%         dot_x = %your code here
%         dot_y = %your code here
%         set(leg_drawing.vertices{vertex_index},'xdata',dot_x,'ydata',dot_y);
%     end
% 
%     %your code here
%     %crank_x and crank_y should both be two element arrays
%     %containing the x and y coordinates of the line segment describing the crank
%     crank_x = %your code here
%     crank_y = %your code here
%     set(leg_drawing.crank,'xdata',crank_x,'ydata',crank_y);
% end