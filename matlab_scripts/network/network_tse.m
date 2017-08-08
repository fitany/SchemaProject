function network = network_tse(flavor_ind,flavor_in_order,network)
    %Schemas (place,flavor)
    schema = [flavor_ind;flavor_in_order]';

    %Neuron sizes
    params.size_wells = 6433; %length(x);%25;
    params.size_flavors = 10;
    params.size_pairs = length(flavor_ind);%5;
    params.size_pfc = 10;
    params.size_hipp = 49;
    params.size_multimodal = 50;
    
    %Learning parameters
    params.lr_slow = .001; %.001
    params.lr_fast = .2; %.2
    params.alpha_multimodal = .2;
    params.alpha_pfc = 5; %2
    params.alpha_hipp = 1; %1
    params.multimodal_excitation = 1;
    params.multimodal_inhibition = 0;
    params.pfc_excitation = 1;
    params.pfc_inhibition = 0;
    params.hipp_excitation = 2;
    params.hipp_inhibition = 0;
    params.gain = 2;
    
    %Training time
    t_per_pair = 3;
    
    %Testing conditions
    has_hipp = 1;
    disp_on = 1;
    sigm = 1;
    
    %Initializing neurons and weights
    if isempty(network)
        network.n_well = zeros(params.size_wells,1); % input neurons for location, 1st half of input
        network.n_flavor = zeros(params.size_flavors,1); % input neurons for flavor, 2nd half of input
        network.n_multimodal = zeros(params.size_multimodal,1);
        network.n_buffer = zeros(numel(network.n_multimodal),1);
        network.n_pfc = rand(params.size_pfc,1).*.1; % output neurons for schema
        network.n_hipp = zeros(params.size_hipp,1); % output neurons for schema
        network.w_input_multimodal = rand(numel(network.n_multimodal),numel(network.n_well)+numel(network.n_flavor)).*.1;
        network.w_buffer_pfc = rand(params.size_pfc,numel(network.n_multimodal)); % weights
        network.w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
        network.w_buffer_hipp = rand(params.size_hipp,numel(network.n_multimodal)).*.1;
    end
    
    %Figure 2A
    network.n_well = zeros(params.size_wells,1); % input neurons for location
    network.n_flavor = zeros(params.size_flavors,1); % input neurons for flavor
    network.n_well(schema(end,1)) = network.n_well(schema(end,1)) + .1;
    network.n_flavor(schema(end,2)) = network.n_flavor(schema(end,2)) + .1;
    network = run_network(network,params,true,t_per_pair,has_hipp,sigm,disp_on);

end