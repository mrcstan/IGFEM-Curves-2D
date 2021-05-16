% test_cosine_min_angle_2D_triangle
coord = [0,0;
         0,1;
         1,1]';
cosMinAngle = cosine_min_angle_2D_triangle(coord)

coord = [0,0;
         0,1;
         0,2]';
cosMinAngle = cosine_min_angle_2D_triangle(coord)

coord = [0,0;
         0,1;
         0.5,sqrt(3)/2.0]';
cosMinAngle = cosine_min_angle_2D_triangle(coord)

coord = [0,0;
         1,0;
         1,0.5]';
cosMinAngle = cosine_min_angle_2D_triangle(coord)