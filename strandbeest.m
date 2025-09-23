function strandbeest()

%Set Parameters
%initialize leg_params structure
leg_params = struct();
%number of vertices in linkage
leg_params.num_vertices = 7;
%number of links in linkage
leg_params.num_linkages = 10;
%matrix relating links to vertices
leg_params.link_to_vertex_list = ...
        [ 1, 3;... %link 1 adjacency
        3, 4;... %link 2 adjacency
        2, 3;... %link 3 adjacency
        2, 4;... %link 4 adjacency
        4, 5;... %link 5 adjacency
        2, 6;... %link 6 adjacency
        1, 6;... %link 7 adjacency
        5, 6;... %link 8 adjacency
        5, 7;... %link 9 adjacency
        6, 7 ... %link 10 adjacency
        ];
leg_params.link_lengths = ...
        [ 50.0,... %link 1 length
        55.8,... %link 2 length
        41.5,... %link 3 length
        40.1,... %link 4 length
        39.4,... %link 5 length
        39.3,... %link 6 length
        61.9,... %link 7 length
        36.7,... %link 8 length
        65.7,... %link 9 length
        49.0 ... %link 10 length
        ];
%length of crank shaft
leg_params.crank_length = 15.0;
%fixed position coords of vertex 0
leg_params.vertex_pos0 = [0;0];
%fixed position coords of vertex 2
leg_params.vertex_pos2 = [-38.0;-7.8];


num_vertices = leg_params.num_vertices;

vertex_coords = %WRITE THE THING HERE;

length_errors = link_length_error_func(vertex_coords, leg_params);

end

function coords_out = column_to_matrix(vertex_coords)

num_coords = length(vertex_coords);
coords_out = [vertex_coords(1:2:(num_coords-1)),vertex_coords(2:2:num_coords)];

end

%Error function that encodes the link length constraints
%INPUTS:
%vertex_coords: a column vector containing the (x,y) coordinates of every vertex
% in the linkage.

%leg_params: a struct containing the parameters that describe the linkage
% importantly, leg_params.link_lengths is a list of linakge lengths
% and leg_params.link_to_vertex_list is a two column matrix where
% leg_params.link_to_vertex_list(i,1) and
% leg_params.link_to_vertex_list(i,2) are the pair of vertices connected
% by the ith link in the mechanism

%OUTPUTS:
%length_errors: a column vector describing the current distance error of the ith
% link specifically, length_errors(i) = (xb-xa)ˆ2 + (yb-ya)ˆ2 - d_iˆ2
% where (xa,ya) and (xb,yb) are the coordinates of the vertices that
% are connected by the ith link, and d_i is the length of the ith link

function length_errors = link_length_error_func(vertex_coords, leg_params)

    link_to_vertex_list = leg_params.link_to_vertex_list;
    link_lengths = leg_params.link_lengths;

    length_errors = zeros(length(leg_params.num_linkages),4);

    for i = 1:leg_params.num_linkages
        vertex_call = link_to_vertex_list(i,:);
        A = vertex_coords(vertex_call(1));
        B = vertex_coords(vertex_call(2));
        fi = (A(1)-B(1))^2 - (A(2)-B(2))^2 - link_lengths(i);

        length_errors(i,:) = fi; 
    end

    disp(length_errors)
end