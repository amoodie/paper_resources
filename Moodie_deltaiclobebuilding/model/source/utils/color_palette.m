function [CP] = color_palette()
%COLOR_PALETTE returns useful pretty 4 color - 5 shade color palette
%
%   c1 = reds
%   c2 = oranges
%   c3 = blues
%   c4 = greens
%
%   s1 = middle, bold shade
%   s2 = lightest
%   s3 = second lightest, lighter than s1
%   s4 = second darkest, darker than s1
%   s5 = darkest
%

    CP.c1_s1 = [0.941,0.408,0.408];
    CP.c1_s2 = [1,0.741,0.741];
    CP.c1_s3 = [1,0.58,0.58];
    CP.c1_s4 = [0.816,0.251,0.251];
    CP.c1_s5 = [0.678,0.129,0.129];
    
    CP.c2_s1 = [0.941,0.651,0.408];
    CP.c2_s2 = [1,0.859,0.741];
    CP.c2_s3 = [1,0.769,0.58];
    CP.c2_s4 = [0.816,0.51,0.251];
    CP.c2_s5 = [0.678,0.376,0.129];
    
    CP.c3_s1 = [0.243,0.565,0.565];
    CP.c3_s2 = [0.62,0.835,0.835];
    CP.c3_s3 = [0.4,0.69,0.69];
    CP.c3_s4 = [0.153,0.49,0.49];
    CP.c3_s5 = [0.078,0.408,0.408];
    
    CP.c4_s1 = [0.325,0.753,0.325];
    CP.c4_s2 = [0.682,0.918,0.682];
    CP.c4_s3 = [0.49,0.847,0.49];
    CP.c4_s4 = [0.2,0.655,0.2];
    CP.c4_s5 = [0.102,0.541,0.102];

    % define some grays
    CP.g1 = 0.2 * ones(1, 3);
    CP.g2 = 0.4 * ones(1, 3);
    CP.g3 = 0.7 * ones(1, 3);
    CP.g4 = 0.9 * ones(1, 3);
    
end