function [stationStruct] = do_velProfCalcs(stationStruct)

    stationStruct.Velocity.Ubar_conCf = sqrt(stationStruct.Velocity.taubDSP / 1000 / stationStruct.Velocity.Cf);
    stationStruct.Velocity.Fr = stationStruct.Velocity.Ubar_conCf / sqrt(9.81 * stationStruct.Velocity.flowDepth);
    
    skin_corr = 'calib';
    if strcmp(skin_corr, 'DSP')
        stationStruct.Velocity.tstar = (stationStruct.Velocity.flowDepth * stationStruct.Velocity.slope) / (1.65 * stationStruct.bedData.gsSummBed.d50.*1e-6);
        stationStruct.Velocity.tstarsk = 0.05 + 0.7 .* (stationStruct.Velocity.tstar .* stationStruct.Velocity.Fr .^ 0.7) .^ 0.8;
        stationStruct.Velocity.ustarsk = sqrt(stationStruct.Velocity.tstarsk * 9.81 * stationStruct.bedData.gsSummBed.d50.*1e-6 .* (1/1.65));
    elseif strcmp(skin_corr, 'calib')
        stationStruct.Velocity.tstar = (stationStruct.Velocity.ustarCalib * stationStruct.Velocity.ustarCalib) / ...
            (9.81 * 1.65 * stationStruct.bedData.gsSummBed.d50.*1e-6);
        stationStruct.Velocity.tstarsk = 0.05 + 0.7 .* (stationStruct.Velocity.tstar .* stationStruct.Velocity.Fr .^ 0.7) .^ 0.8;
        stationStruct.Velocity.ustarsk = sqrt(stationStruct.Velocity.tstarsk * 9.81 * stationStruct.bedData.gsSummBed.d50.*1e-6 .* (1.65));
    end

    
    % Ma friction relation
    if ~isempty(stationStruct.Velocity.data.measZ)
        Usq = (stationStruct.Velocity.data.mean).^2;
        ws = get_DSV(stationStruct.bedData.gsSummBed.d50 ./ 1e6, 0.7, 3.5, load_conset('quartz-water'));
        dustar = 1;
        ustar_fric = stationStruct.Velocity.ustarDSP;
        iter = 0;
        while dustar > 0.000001
            usws = ustar_fric / ws;
            y = -1.53 * (log10(usws)^2) - 0.15 * (log10(usws)) - 1.9;
            newus = sqrt( Usq * 10^(y) );
            differ = newus - ustar_fric;
            dustar = differ / 8;
            ustar_fric = ustar_fric + dustar;
            iter = iter + 1;
        end
        stationStruct.Velocity.ustarFric = ustar_fric;
        stationStruct.Velocity.ustarFric_MaSolved = Cf_Ma_getUs(stationStruct.Velocity.data.mean, ws);
    end
end
