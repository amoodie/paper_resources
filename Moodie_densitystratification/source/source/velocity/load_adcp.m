function [adcp] = load_adcp(stnInfoMod, tddepth)

    data_filename = strcat('./datasrc/adcp_ascii/', stnInfoMod, '_data.TXT');
    depth_filename = strcat('./datasrc/adcp_ascii/', stnInfoMod, '_depth.TXT');
    
    adcp.data = csvread(data_filename, 0, 0);
    
    adcp.data(adcp.data == -32768) = NaN; % replace bad data with NaN
    
    adcp.depth = csvread(depth_filename, 0, 0);
    
    adcp.tddepth = tddepth; % transducer depth
    
end