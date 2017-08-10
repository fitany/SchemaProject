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
    params.lr_slow = .01; %.01 %.001
    params.lr_fast = .3; %.3 %.2
    params.alpha_multimodal = .2;
    params.alpha_pfc = 2; %2
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
    num_trials = 20;
    performances = zeros(1,num_trials);
    for trial = 1:num_trials
        network = train_schema(schemaA,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
        performances(trial) = test_schema(schemaA,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    end
    %for trial = 1:num_trials
    %    network = train_schema(schemaB,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    %end
    %load('trained_network.mat');
    %performance = test_schema(schemaA,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    figure;
    plot(1:num_trials,performances,'o-');
    title('Performance over time training on schema A');
    xlabel('Day');
    ylabel('Performance');
end
function performance = test_schema(schema,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
    performance = 0;
    %Run network without plasticity to find the active schema
    for p = 1:params.size_pairs;
    %for p = 1:t_pairs;
        network.n_well = zeros(params.size_wells,1); % input neurons for location
        network.n_flavor = zeros(params.size_flavors,1); % input neurons for flavor
        %i = randi(params.size_pairs);
        i=p;
        network.n_well(schema(i,1)) = network.n_well(schema(i,1)) + .1;
        network.n_flavor(schema(i,2)) = network.n_flavor(schema(i,2)) + .1;
        network = run_network(network,params,false,t_per_pair,has_hipp,sigm,disp_on);
    end
    %Find active schema and set pfc to that
    [~,mInd] = max(network.n_pfc);
    network.n_pfc = network.n_pfc .* 0;
    network.n_pfc(mInd) = 1;
    %Cue each flavor to find associated location
    for p = 1:params.size_pairs
        network = reset(network,params,true,false);
        network.n_pfc = zeros(params.size_pfc,1);
        network.n_pfc(mInd) = 1;
        network.n_well = zeros(params.size_wells,1); % input neurons for location
        network.n_flavor = zeros(params.size_flavors,1); % input neurons for flavor
        i = p;
        network.n_flavor(schema(i,2)) = network.n_flavor(schema(i,2)) + 1;
        network = run_network_retrieval_simple(network,params,1,has_hipp,sigm,disp_on);
        well_locations = schema(:,1);
        all_predictions = network.n_well./sum(network.n_well);
        well_predictions = all_predictions(well_locations)./sum(all_predictions(well_locations));
        well_predictions = exp(100.*well_predictions)./sum(exp(100.*well_predictions)); % applied softmax
        performance = performance + well_predictions(i);
        fprintf('Expected location: %d, Predictions: %d:%1.2f,%d:%1.2f,%d:%1.2f,%d:%1.2f,%d:%1.2f\n',schema(i,1),well_locations(1),well_predictions(1),well_locations(2),well_predictions(2),well_locations(3),well_predictions(3),well_locations(4),well_predictions(4),well_locations(5),well_predictions(5))
    end
    performance = performance./params.size_pairs
end

function network = train_schema(schema,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_flavors = params.size_flavors;
    pair_order = randperm(params.size_pairs);
    for p = 1:params.size_pairs
        network.n_well = zeros(size_wells,1); % input neurons for location
        network.n_flavor = zeros(size_flavors,1); % input neurons for flavor
        i=pair_order(p);
        network.n_well(schema(i,1)) = network.n_well(schema(i,1)) + .1;
        network.n_flavor(schema(i,2)) = network.n_flavor(schema(i,2)) + .1;
        network = run_network(network,params,true,t_per_pair,has_hipp,sigm,disp_on);
    end
end

function network = reset(network,params,reset_pop,reset_weights)
    if reset_pop
        network.n_well = zeros(params.size_wells,1); % input current for location, 1st half of input
        network.n_flavor = zeros(params.size_flavors,1); % input current for flavor, 2nd half of input
        network.n_input = zeros(numel(network.n_well)+numel(network.n_flavor),1); % layer receiving input
        network.n_multimodal = zeros(params.size_multimodal,1);
        network.n_buffer = zeros(numel(network.n_multimodal),1);
        network.n_pfc = rand(params.size_pfc,1).*.1; % output neurons for schema
        network.n_hipp = zeros(params.size_hipp,1); % output neurons for schema
    end
    if reset_weights
        network.w_input_multimodal = rand(numel(network.n_multimodal),numel(network.n_well)+numel(network.n_flavor)).*.1;
        network.w_buffer_pfc = rand(params.size_pfc,numel(network.n_multimodal)); % weights
        network.w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
        network.w_buffer_hipp = rand(params.size_hipp,numel(network.n_multimodal)).*.1;
    end
end
