% Old code for reference
% Final Project
% Tiffany Hwu
% ------------------------------------------------------------------------------------
function tse_hipp_wta_DA_buffer()
    %Schemas
    schema1 = [9,1;14,2;19,3;31,4;36,5;41,6];
    schema2 = [3,11;20,12;24,13;29,14;33,15;47,16];
    schema3 = [3,15;20,14;24,12;29,16;33,11;47,13];
    
    %Neuron sizes
    params.size_wells = 49;
    params.size_flavors = 16;
    params.size_pairs = 6;
    params.size_pfc = 10;
    params.size_hipp = 49;
    
    %Learning parameters
    params.lr_slow = .001; %.001
    params.lr_fast = .2; %.2
    params.alpha_pfc = 5; %2
    params.alpha_hipp = 1; %1
    params.pfc_excitation = 1;
    params.pfc_inhibition = 0;
    params.hipp_excitation = 2;
    params.hipp_inhibition = 0;
    params.gain = 2;
    
    %Training time
    t_per_pair = 10;
    t_pairs = 50;
    
    %Testing conditions
    has_hipp = 1;
    disp_on = 1;
    sigm = 0;
    
    %Initializing neurons and weights
    network.n_well = zeros(params.size_wells,params.size_pairs); % input neurons for location
    network.n_flavor = zeros(params.size_flavors,params.size_pairs); % input neurons for flavor
    network.n_input = zeros(numel(network.n_well)+numel(network.n_flavor),1);
    network.n_buffer = zeros(numel(network.n_well)+numel(network.n_flavor),1);
    network.n_pfc = rand(params.size_pfc,1); % output neurons for schema
    network.n_hipp = zeros(params.size_hipp,1); % output neurons for schema
    network.w_sb_pfc = rand(params.size_pfc,numel(network.n_well)+numel(network.n_flavor)); % weights
    network.w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
    network.w_sc_hipp = zeros(params.size_hipp,numel(network.n_well)+numel(network.n_flavor));
    
    %[w_sb_pfc,w_pfc_hipp,w_sc_hipp,n_input,n_buffer,n_pfc,n_hipp] = train_schema(schema1,n_input,n_buffer,n_pfc,n_hipp,w_sb_pfc,w_pfc_hipp,w_sc_hipp,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    %n_input = zeros(numel(n_well)+numel(n_flavor),1);
    %n_buffer = zeros(numel(n_well)+numel(n_flavor),1);
    %n_hipp = zeros(params.size_hipp,1);
    %w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
    %w_sc_hipp = zeros(params.size_hipp,numel(n_well)+numel(n_flavor));
    %performance = test_schema(schema1,n_input,n_buffer,n_pfc,n_hipp,w_sb_pfc,w_pfc_hipp,w_sc_hipp,params,t_per_pair,has_hipp,sigm,disp_on)
    
    %disp 'Training Schema 2'
    %n_input = zeros(numel(n_well)+numel(n_flavor),1);
    %n_buffer = zeros(numel(n_well)+numel(n_flavor),1);
    %n_hipp = zeros(params.size_hipp,1);
    %w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
    %w_sc_hipp = zeros(params.size_hipp,numel(n_well)+numel(n_flavor));
    %[w_sb_pfc,w_pfc_hipp,w_sc_hipp,n_input,n_buffer,n_pfc,n_hipp] = train_schema(schema2,n_input,n_buffer,n_pfc,n_hipp,w_sb_pfc,w_pfc_hipp,w_sc_hipp,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
    
    %Figure 2A
    num_trials = 19;
    performances = zeros(1,num_trials);
    for trial = 1:num_trials
        network = train_schema(schema1,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on);
        test_network.n_input = zeros(numel(network.n_well)+numel(network.n_flavor),1);
        test_network.n_buffer = zeros(numel(network.n_well)+numel(network.n_flavor),1);
        test_network.n_hipp = zeros(params.size_hipp,1);
        test_network.w_pfc_hipp = zeros(params.size_hipp,params.size_pfc);
        test_network.w_sc_hipp = zeros(params.size_hipp,numel(network.n_well)+numel(network.n_flavor));
        test_network.w_sb_pfc = network.w_sb_pfc;
        test_network.n_pfc = network.n_pfc;
        performance = test_schema(schema1,test_network,params,t_per_pair,has_hipp,sigm,disp_on)
        performances(trial) = performance;
    end
    performances
end

function network = train_schema(schema,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_flavors = params.size_flavors;
    size_pairs = params.size_pairs;
    size_pfc = params.size_pfc;
    size_hipp = params.size_hipp;
    
    lr_slow = params.lr_slow;
    lr_fast = params.lr_fast;
    alpha_pfc = params.alpha_pfc;
    alpha_hipp = params.alpha_hipp;
    pfc_excitation = params.pfc_excitation;
    pfc_inhibition = params.pfc_inhibition;
    hipp_excitation = params.hipp_excitation;
    hipp_inhibition = params.hipp_inhibition;
    gain = params.gain;
    
    w_sb_pfc = network.w_sb_pfc;
    w_pfc_hipp = network.w_pfc_hipp;
    w_sc_hipp = network.w_sc_hipp;
    n_input = network.n_input;
    n_buffer = network.n_buffer;
    n_pfc = network.n_pfc;
    n_hipp = network.n_hipp;
    
    pair_order = randperm(size_pairs);
    
    for p = 1:size_pairs
    %for p = 1:t_pairs
        %Input current for schema 1
        n_well = zeros(size_wells,size_pairs); % input neurons for location
        n_flavor = zeros(size_flavors,size_pairs); % input neurons for flavor
        %i = randi(params.size_pairs);
        i = pair_order(p);
        n_well(schema(i,1),i) = n_well(schema(i,1),i) + 1;
        n_flavor(schema(i,2),i) = n_flavor(schema(i,2),i) + 1;
        input_current = [reshape(n_well,numel(n_well),1);reshape(n_flavor,numel(n_flavor),1)];
        
        for t = 1:t_per_pair
            %Update input neurons
            if has_hipp
                if sigm
                    n_input_new = sigmoid(w_sc_hipp' * n_hipp + input_current,gain1)-.5; %activation
                else
                    n_input_new = max(w_sc_hipp' * n_hipp + input_current,0); %activation
                end
            else
                if sigm
                    n_input_new = sigmoid(input_current,gain)-.5; %activation
                else
                    n_input_new = max(input_current,0); %activation
                end
            end
            %Update input buffer neurons
            n_buffer = .99 .* n_buffer;
            n_buffer = max([n_input_new'; n_buffer'; zeros(size(n_input))'])';
            %Update pfc neurons
            if has_hipp
                if sigm
                    n_pfc_new = sigmoid(w_sb_pfc * n_buffer + w_pfc_hipp' * n_hipp,gain)-.5; %activation
                else
                    n_pfc_new = max(w_sb_pfc * n_buffer + w_pfc_hipp' * n_hipp,0); %activation
                end
            else
                n_pfc_new = sigmoid(w_sb_pfc * n_buffer,gain)-.5; %activation
            end
            [~,mInd] = max(n_pfc_new);
            pfc_activity =  n_pfc_new(mInd);
            %if pfc_activity > 50
            %    disp 'schema active';
            %    lr_fast = .4;
            %    lr_slow = .1;
            %end
            n_pfc_new = -pfc_inhibition * ones(size_pfc,1);
            n_pfc_new(mInd) = pfc_excitation;
            %Update hipp neurons
            if has_hipp
                if sigm
                    n_hipp_new = sigmoid(w_pfc_hipp * n_pfc + w_sc_hipp * n_input,gain)-.5; %activation
                else
                    n_hipp_new = max(w_pfc_hipp * n_pfc + w_sc_hipp * n_input,0); %activation
                end
                mInd = schema(i,1); % external gps signal
                n_hipp_new = -hipp_inhibition * ones(size_hipp,1);
                n_hipp_new(mInd) = hipp_excitation;
            end

            %Update weights
            w_sb_pfc = w_sb_pfc + ojaLearningRule(lr_slow, n_buffer, n_pfc, w_sb_pfc, alpha_pfc);
            if has_hipp
                w_pfc_hipp = w_pfc_hipp + ojaLearningRule(lr_fast, n_pfc, n_hipp, w_pfc_hipp, alpha_hipp);
                w_sc_hipp = w_sc_hipp + ojaLearningRuleDA(lr_fast, n_input, n_hipp, w_sc_hipp, alpha_hipp, find(input_current));
            end

            if disp_on
                subplot(2,4,1);
                imagesc(n_input);
                title('n input');
                subplot(2,4,2)
                imagesc(n_buffer);
                title('n buffer');
                subplot(2,4,3)
                imagesc(n_pfc);
                title('n pfc');
                subplot(2,4,4)
                imagesc(n_hipp);
                title('n hipp');       
                subplot(2,4,5);
                imagesc(w_sb_pfc);
                title('w s pfc');
                subplot(2,4,6);
                imagesc(w_pfc_hipp);
                title('w pfc hipp');
                subplot(2,4,7);
                imagesc(w_sc_hipp);
                title('w s hipp');
                pause(.001);
            end

            n_input = n_input_new;
            n_pfc = n_pfc_new;
            if has_hipp
                n_hipp = n_hipp_new;
            end
        end
    end
    % test effectiveness of training
    % sum(w_sb_pfc(:,[63,302,117,309,178,316,232,323,286,330,9,295]),2)./sum(w_sb_pfc,2)
    
    network.w_sb_pfc = w_sb_pfc;
    network.w_pfc_hipp = w_pfc_hipp;
    network.w_sc_hipp = w_sc_hipp;
    network.n_input = n_input;
    network.n_buffer = n_buffer;
    network.n_pfc = n_pfc;
    network.n_hipp = n_hipp;
end

function network = train_new_PA(schema,removed_PAs,new_PAs,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
    %Remove PAs
    for i = 1:size(removed_PAs,1)
        for j = 1:size(schema,1)
            if sum(schema(j,:) == removed_PAs(i,:)) == 2
                schema(j,:) = [];
            end
        end
    end
    %Add PAs
    schema = [schema; new_PAs];

    network = train_schema(schema,network,params,t_per_pair,t_pairs,has_hipp,sigm,disp_on)
end

%add performance calculation
function performance = test_schema(schema,network,params,t_per_pair,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_flavors = params.size_flavors;
    size_pairs = params.size_pairs;
    size_pfc = params.size_pfc;
    size_hipp = params.size_hipp;
    
    lr_slow = params.lr_slow;
    lr_fast = params.lr_fast;
    alpha_pfc = params.alpha_pfc;
    alpha_hipp = params.alpha_hipp;
    pfc_excitation = params.pfc_excitation;
    pfc_inhibition = params.pfc_inhibition;
    hipp_excitation = params.hipp_excitation;
    hipp_inhibition = params.hipp_inhibition;
    gain = params.gain;
    
    w_sb_pfc = network.w_sb_pfc;
    w_pfc_hipp = network.w_pfc_hipp;
    w_sc_hipp = network.w_sc_hipp;
    n_input = network.n_input;
    n_buffer = network.n_buffer;
    n_pfc = network.n_pfc;
    n_hipp = network.n_hipp;
    
    performances = zeros(1,size(schema,1));
    
    for p = 1:size(schema,1)
        %Input current for schema 1
        n_well = zeros(params.size_wells,params.size_pairs); % input neurons for location
        n_flavor = zeros(params.size_flavors,params.size_pairs); % input neurons for flavor
        i = p;
        n_flavor(schema(i,2),i) = n_flavor(schema(i,2),i) + 1;
        input_current = [reshape(n_well,numel(n_well),1);reshape(n_flavor,numel(n_flavor),1)];
        
        for t = 1:t_per_pair
            %Update input neurons
            if has_hipp
                if sigm
                    n_input_new = sigmoid(w_sc_hipp' * n_hipp + input_current + w_sb_pfc' * n_pfc,gain); %activation
                else
                    n_input_new = max(w_sc_hipp' * n_hipp + input_current + w_sb_pfc' * n_pfc,0); %activation
                end
            else
                if sigm
                    n_input_new = sigmoid(input_current + w_sb_pfc' * n_pfc,gain); %activation
                else
                    n_input_new = max(input_current + w_sb_pfc' * n_pfc,0); %activation
                end
            end
            %Update input buffer neurons
            n_buffer = .99 .* n_buffer;
            n_buffer = max([n_input_new'; n_buffer'; zeros(size(n_input))'])';
            %Update pfc neurons
            if has_hipp
                if sigm
                    n_pfc_new = sigmoid(w_sb_pfc * n_buffer + w_pfc_hipp' * n_hipp,gain); %activation
                else
                    n_pfc_new = max(w_sb_pfc * n_buffer + w_pfc_hipp' * n_hipp,0); %activation
                end
            else
                if sigm
                    n_pfc_new = sigmoid(w_sb_pfc * n_buffer,gain); %activation
                else
                    n_pfc_new = max(w_sb_pfc * n_buffer,0); %activation
                end
            end
            [~,mInd] = max(n_pfc_new);
            %pfc_excitation = n_pfc_new(mInd);
            n_pfc_new = -pfc_inhibition * ones(size_pfc,1);
            n_pfc_new(mInd) = pfc_excitation;
            %Update hipp neurons
            if has_hipp
                if sigm
                    n_hipp_new = sigmoid(w_pfc_hipp * n_pfc + w_sc_hipp * n_input,gain); %activation
                else
                    n_hipp_new = max(w_pfc_hipp * n_pfc + w_sc_hipp * n_input,0); %activation
                end
                mInd = schema(i,1); % external gps signal
                n_hipp_new = -hipp_inhibition * ones(size_hipp,1);
                n_hipp_new(mInd) = hipp_excitation;
            end

            if disp_on
                subplot(2,4,1);
                imagesc(n_input);
                title('n input');
                subplot(2,4,2)
                imagesc(n_buffer);
                title('n buffer');
                subplot(2,4,3)
                imagesc(n_pfc);
                title('n pfc');
                subplot(2,4,4)
                imagesc(n_hipp);
                title('n hipp');       
                subplot(2,4,5);
                imagesc(w_sb_pfc);
                title('w s pfc');
                subplot(2,4,6);
                imagesc(w_pfc_hipp);
                title('w pfc hipp');
                subplot(2,4,7);
                imagesc(w_sc_hipp);
                title('w s hipp');
                pause(.001);
            end

            n_input = n_input_new;
            n_pfc = n_pfc_new;
            if has_hipp
                n_hipp = n_hipp_new;
            end
        end
        %Reshaping input back into well and flavor
        n_well = reshape(n_input(1:numel(n_well)),size_wells,size_pairs);
        n_flavor = reshape(n_input(numel(n_well)+1:end),size_flavors,size_pairs);
        wells = n_well(schema(:,1),i);
        performances(p) = wells(p)/sum(wells);
    end
    performance = mean(performances);
end

function [performanceCue,performanceNonCue] = probe_test(new_PAs,n_input,n_buffer,n_pfc,n_hipp,w_sb_pfc,w_pfc_hipp,w_sc_hipp,params,t_per_pair,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_flavors = params.size_flavors;
    size_pairs = params.size_pairs;
    size_pfc = params.size_pfc;
    size_hipp = params.size_hipp;
    
    lr_slow = params.lr_slow;
    lr_fast = params.lr_fast;
    alpha_pfc = params.alpha_pfc;
    alpha_hipp = params.alpha_hipp;
    pfc_excitation = params.pfc_excitation;
    pfc_inhibition = params.pfc_inhibition;
    hipp_excitation = params.hipp_excitation;
    hipp_inhibition = params.hipp_inhibition;
    gain = params.gain;
    
    %Input current for PA
    n_well = zeros(params.size_wells,params.size_pairs); % input neurons for location
    n_flavor = zeros(params.size_flavors,params.size_pairs); % input neurons for flavor
    for i = 1:size(schema,1)
        if sum(schema(i,:) == new_PAs(1,:)) == 2
            break;
        end
    end
    for j = 1:size(schema,1)
        if sum(schema(i,:) == new_PAs(2,:)) == 2
            break;
        end
    end
    n_flavor(schema(i,2),i) = n_flavor(schema(i,2),i) + 1;
    input_current = [reshape(n_well,numel(n_well),1);reshape(n_flavor,numel(n_flavor),1)];
        
    for t = 1:t_per_pair
        %Update input neurons
        if has_hipp
            if sigm
                n_input_new = sigmoid(w_sc_hipp' * n_hipp + input_current,gain); %activation
            else
                n_input_new = max(w_sc_hipp' * n_hipp + input_current,0); %activation
            end
        else
            if sigm
                n_input_new = sigmoid(input_current,gain); %activation
            else
                n_input_new = max(input_current,0); %activation
            end
        end
        %Update input buffer neurons
        n_buffer = .99 .* n_buffer;
        n_buffer = max([n_input_new'; n_buffer'; zeros(size(n_input))'])';
        %Update pfc neurons
        if has_hipp
            if sigm
                n_pfc_new = sigmoid(w_sb_pfc * n_buffer + w_pfc_hipp' * n_hipp,gain); %activation
            else
                n_pfc_new = max(w_sb_pfc * n_buffer + w_pfc_hipp' * n_hipp,0); %activation
            end
        else
            if sigm
                n_pfc_new = sigmoid(w_sb_pfc * n_buffer,gain); %activation
            else
                n_pfc_new = max(w_sb_pfc * n_buffer,0); %activation
            end
        end
        [~,mInd] = max(n_pfc_new);
        n_pfc_new = -pfc_inhibition * ones(size_pfc,1);
        n_pfc_new(mInd) = pfc_excitation;
        %Update hipp neurons
        if has_hipp
            if sigm
                n_hipp_new = sigmoid(w_pfc_hipp * n_pfc + w_sc_hipp * n_input,gain); %activation
            else
                n_hipp_new = max(w_pfc_hipp * n_pfc + w_sc_hipp * n_input,0); %activation
            end
            mInd = schema(i,1); % external gps signal
            n_hipp_new = -hipp_inhibition * ones(size_hipp,1);
            n_hipp_new(mInd) = hipp_excitation;
        end

        if disp_on
            subplot(2,4,1);
            imagesc(n_input);
            title('n input');
            subplot(2,4,2)
            imagesc(n_buffer);
            title('n buffer');
            subplot(2,4,3)
            imagesc(n_pfc);
            title('n pfc');
            subplot(2,4,4)
            imagesc(n_hipp);
            title('n hipp');       
            subplot(2,4,5);
            imagesc(w_sb_pfc);
            title('w s pfc');
            subplot(2,4,6);
            imagesc(w_pfc_hipp);
            title('w pfc hipp');
            subplot(2,4,7);
            imagesc(w_sc_hipp);
            title('w s hipp');
            pause(.001);
        end

        n_input = n_input_new;
        n_pfc = n_pfc_new;
        if has_hipp
            n_hipp = n_hipp_new;
        end
    end
    
    %Reshaping input back into well and flavor
    n_well = reshape(n_input(1:numel(n_well)),size_wells,size_pairs);
    n_flavor = reshape(n_input(numel(n_well)+1:end),size_flavors,size_pairs);
    wells = n_well(schema(:,1),i);
    performance_cue = wells(i)/sum(wells);
    performance_non_cue = wells(j)/sum(wells);
end

function gain = sigmoid(x,k)
    gain = 1./(1+exp(-k*x));
end

function dW = ojaLearningRule(lr, pre, post, w, alpha)
    dW = lr .* ((post*pre') - alpha.*repmat(post.^2,1,size(w,2)).*w);
end

function dW = ojaLearningRuleDA(lr, pre, post, w, alpha, DA)
    reward = zeros(size(pre));
    reward(DA) = 1;
    dW = lr .* ((post*(pre .* reward)') - alpha.*repmat(post.^2,1,size(w,2)).*w);
end