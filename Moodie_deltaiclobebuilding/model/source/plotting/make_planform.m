function [plan] = make_planform(x, gamm, gammapex, rad, addapex)
    plan = zeros(100, 2);
    pts = linspace( deg2rad(90), deg2rad(90-gamm), 100);
    plan(:, 1) = cosd(gamm/2)*gammapex + ((rad-gammapex) * cos(pts)); % x coords
    plan(:, 2) = sind(gamm/2)*gammapex + ((rad-gammapex) * sin(pts)); % y coords
    
    if addapex
        plan(end+1, :) = [cosd(gamm/2)*gammapex, sind(gamm/2)*gammapex];
    end
    
end
