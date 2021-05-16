% test_cosine_min_angle_2D_triangle
%{
coord = [0,0;
         0,1;
         1,1]';
ratio = aspect_ratio_2D(coord)

coord = [0,0;
         0,1;
         0,2]';
ratio = aspect_ratio_2D(coord)

coord = [0,0;
         0,1;
         0.5,sqrt(3)/2.0]';
ratio = aspect_ratio_2D(coord)

coord = [0,0;
         1,0;
         1,0.5]';
ratio = aspect_ratio_2D(coord)

coord = [0,0;
         1,0;
         1,1;
         0,1]';
ratio = aspect_ratio_2D(coord)   
%}
coord = [0,0;
         1,0;
         1,1;
         0,0.01]';
ratio = aspect_ratio_2D(coord)  