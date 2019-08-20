function [river, lobe] = make_lobeplan(x, gamm, rad, mou_idx, Bo0, dx)
    mouxy = [cosd(gamm/2)*mou_idx*dx, sind(gamm/2)*mou_idx*dx];
    radxy = [cosd(gamm/2)*(rad-(1*dx)), sind(gamm/2)*(rad-(1*dx))];
    river = [0, 0; mouxy];

    lobe = [radxy(1) + sind(gamm/2)*(Bo0/2), radxy(2) - cosd(gamm/2)*(Bo0/2);
            mouxy(1) + sind(gamm/2)*(Bo0/2), mouxy(2) - cosd(gamm/2)*(Bo0/2);
            mouxy(1) - sind(gamm/2)*(Bo0/2), mouxy(2) + cosd(gamm/2)*(Bo0/2);
            radxy(1) - sind(gamm/2)*(Bo0/2), radxy(2) + cosd(gamm/2)*(Bo0/2);];

end