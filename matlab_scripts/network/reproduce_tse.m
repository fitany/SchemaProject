function reproduce_tse()
    %Schemas (place,flavor)
    schemaA = [4,4;6,1;13,5;15,3;22,2];
    schemaA_swapped = [4,4;6,1;13,5;15,3;22,6];
    schemaA_moved = [4,4;6,1;13,5;20,3;22,2];
    schemaB = [4,7;6,1;8,9;19,8;25,10];
    schemaC = [3,10;10,1;11,2;16,8;24,7];
    
    %Neuron sizes
    params.size_wells = 25;
    params.size_flavors = 10;
    params.size_pairs = 5;
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
    t_pairs = 20;
    
    %Testing conditions
    has_hipp = 1;
    disp_on = 1;
    sigm = 1;
    
    %Initializing neurons and weights
    network.n_well = zeros(params.size_wells,1); % input current for location, 1st half of input
    network.n_flavor = zeros(params.size_flavors,1); % input current for flavor, 2nd half of input
    network.n_input = zeros(numel(network.n_well)+numel(network.n_flavor),1); % layer receiving input
    network.n_multimodal = zeros(params.size_multimodal,1);
    network.n_buffer = zeros(numel(network.n_multimodal),1);
    network.n_pfc = rand(params.size_pfc,1).*.1; % output neurons for schema
    network.n_hipp = zeros(params.size_hipp,1); % output neurons for schema
    network.w_input_multimodal = rand(numel(network.n_multimodal),numel(network.n_well)+numel(network.n_flavor)).*.1;
    network.w_buffer_pfc = rand(params.size_pfc,numel(network.n_multimodal)); % weights
    network.w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
    network.w_buffer_hipp = rand(params.size_hipp,numel(network.n_multimodal)).*.1;
    
    %Figure 2A
    num_trials = 19;%19
    performances = zeros(1,num_trials);
    for trial = 1:num_trials
        network = train_schema(schemaA,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
        %{
        test_network.n_multimodal = zeros(numel(network.n_well)+numel(network.n_flavor),1);
        test_network.n_buffer = zeros(numel(network.n_well)+numel(network.n_flavor),1);
        test_network.n_hipp = zeros(params.size_hipp,1);
        test_network.w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
        test_network.w_buffer_hipp = zeros(params.size_hipp,numel(network.n_well)+numel(network.n_flavor));
        test_network.w_buffer_pfc = network.w_buffer_pfc;
        test_network.n_pfc = network.n_pfc;
        performance = test_schema(schemaA,test_network,params,t_per_pair,has_hipp,sigm,disp_on)
        performances(trial) = performance;
        %}
    end
    %for trial = 1:num_trials
    %    network = train_schema(schemaB,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    %end
    for trial = 1:num_trials
        network = test_schema(schemaA,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    end
    %performances
end
function network = test_schema(schema,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_flavors = params.size_flavors;
    size_pairs = params.size_pairs;
    for p = 1:t_pairs
        n_well = zeros(size_wells,1); % input neurons for location
        n_flavor = zeros(size_flavors,1); % input neurons for flavor
        i = randi(params.size_pairs);
        n_well(schema(i,1)) = n_well(schema(i,1)) + .1;
        n_flavor(schema(i,2)) = n_flavor(schema(i,2)) + .1;
        network.n_well = n_well;
        network.n_flavor = n_flavor;
        network = run_network(n_well,n_flavor,network,params,true,t_per_pair,has_hipp,sigm,disp_on);
    end
    for p = 1:size_pairs
        n_well = zeros(size_wells,1); % input neurons for location
        n_flavor = zeros(size_flavors,1); % input neurons for flavor
        i = p;
        n_flavor(schema(i,2)) = n_flavor(schema(i,2)) + .1;
        network.n_well = n_well;
        network.n_flavor = n_flavor;
        network = run_network_retrieval_simple(network,params,t_per_pair,has_hipp,sigm,disp_on);
    end
end

function network = train_schema(schema,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_flavors = params.size_flavors;
    for p = 1:t_pairs
        n_well = zeros(size_wells,1); % input neurons for location
        n_flavor = zeros(size_flavors,1); % input neurons for flavor
        i = randi(params.size_pairs);
        n_well(schema(i,1)) = n_well(schema(i,1)) + .1;
        n_flavor(schema(i,2)) = n_flavor(schema(i,2)) + .1;
        network.n_well = n_well;
        network.n_flavor = n_flavor;
        network = run_network(n_well,n_flavor,network,params,true,t_per_pair,has_hipp,sigm,disp_on);
    end
end

