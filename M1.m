%M1: run 2_layermodel first to set initial model setup


gravity_anomaly = img2mx
for rho = 2500:10:2800;

    Model.l1.dens = rho
    %%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Part that cant be modified %%%%%%%%%%%%%%%%%%%%%%%

    latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
    lonLim =    [-180 180 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
    height =    25.0; % height of computation above spheroid
    SHbounds =  [0 100]; % Truncation settings: lower limit, upper limit SH-coefficients used

    %%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

    %% Global Spherical Harmonic Analysis 

    tic;
    [V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
    toc

    %% Global Spherical Harmonic Synthesis

    tic;
    [data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_Model,Model);
    toc

    %% Save data

    save([HOME '/Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '.mat'],'data','V_Model','Model')

end