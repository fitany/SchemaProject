% Reverses the weights, input->pfc->buffer->multimodal->well and flavor
% Input with flavor and schema to retrieve location
function network = run_network_retrieval_simple(network,params,t_per_pair,has_hipp,sigm,disp_on)
    size_wells = params.size_wells;
    size_pfc = params.size_pfc;
    size_hipp = params.size_hipp;
    size_multimodal = params.size_multimodal;
    
    lr_slow = params.lr_slow;
    lr_fast = params.lr_fast;
    alpha_multimodal = params.alpha_multimodal;
    alpha_pfc = params.alpha_pfc;
    alpha_hipp = params.alpha_hipp;
    multimodal_excitation = params.multimodal_excitation;
    multimodal_inhibition = params.multimodal_inhibition;
    pfc_excitation = params.pfc_excitation;
    pfc_inhibition = params.pfc_inhibition;
    hipp_excitation = params.hipp_excitation;
    hipp_inhibition = params.hipp_inhibition;
    gain = params.gain;
    
    w_input_multimodal = network.w_input_multimodal;
    w_buffer_pfc = network.w_buffer_pfc;
    w_pfc_hipp = network.w_pfc_hipp;
    w_buffer_hipp = network.w_buffer_hipp;
    n_well = network.n_well;
    n_flavor = network.n_flavor;
    n_multimodal = network.n_multimodal;
    n_buffer = network.n_buffer;
    n_pfc = network.n_pfc;
    n_hipp = network.n_hipp;
    
    %Input current for schema 1
    input_current = [n_well; n_flavor];

    for t = 1:t_per_pair
        %Update input buffer neurons
        if sigm
            n_buffer = 2.*sigmoid(w_buffer_pfc'*n_pfc,gain)-1; %activation
        else
            n_buffer = max(w_buffer_pfc'*n_pfc,0); %activation
        end
        %Update multimodal neurons
        if sigm
            n_multimodal = 2.*sigmoid(w_input_multimodal*input_current + n_buffer,gain)-1; %activation
        else
            n_multimodal = max(w_input_multimodal*input_current + n_buffer,0); %activation
        end
        %Update input neurons
        if sigm
            n_input = 2.*sigmoid(w_input_multimodal'*n_multimodal,gain)-1; %activation
        else
            n_input = max(w_input'*n_multimodal,0); %activation
        end
        
        n_well = n_input(1:size_wells);
        n_flavor = n_input(size_wells+1:end);

        if disp_on
            subplot(3,3,1);
            imagesc([n_well;n_flavor]);
            title('n input');
            subplot(3,3,2);
            imagesc(n_multimodal);
            title('n multimodal');
            subplot(3,3,3)
            imagesc(n_buffer);
            title('n buffer');
            subplot(3,3,4)
            imagesc(n_pfc);
            title('n pfc');
            subplot(3,3,5)
            imagesc(n_hipp);
            title('n hipp'); 
            subplot(3,3,6);
            imagesc(w_input_multimodal);
            title('w input multimodal');
            subplot(3,3,7);
            imagesc(w_buffer_pfc);
            title('w buffer pfc');
            subplot(3,3,8);
            imagesc(w_pfc_hipp);
            title('w pfc hipp');
            subplot(3,3,9);
            imagesc(w_buffer_hipp);
            title('w buffer hipp');
            pause(.001);
        end
    end

    network.w_input_multimodal = w_input_multimodal;
    network.w_buffer_pfc = w_buffer_pfc;
    network.w_pfc_hipp = w_pfc_hipp;
    network.w_buffer_hipp = w_buffer_hipp;
    network.n_well = n_well;
    network.n_flavor = n_flavor;
    network.n_multimodal = n_multimodal;
    network.n_buffer = n_buffer;
    network.n_pfc = n_pfc;
    network.n_hipp = n_hipp;
end

function gain = sigmoid(x,k)
    gain = 1./(1+exp(-k*x));
end

function dW = ojaLearningRule(lr, pre, post, w, alpha)
    dW = lr .* ((post*pre') - alpha.*repmat(post.^2,1,size(w,2)).*w);
end