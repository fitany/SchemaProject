function network = run_network(n_well,n_flavor,network,params,is_learning,t_per_pair,has_hipp,sigm,disp_on)
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
    n_multimodal = network.n_multimodal;
    n_buffer = network.n_buffer;
    n_pfc = network.n_pfc;
    n_hipp = network.n_hipp;
    
    %Input current for schema 1
    input_current = [n_well; n_flavor];

    for t = 1:t_per_pair
        %Update multimodal neurons
        if sigm
            n_multimodal_new = 2.*sigmoid(w_input_multimodal*input_current,gain)-1; %activation
        else
            n_multimodal_new = max(w_input_multimodal*input_current,0); %activation
        end
        [max_multimodal,mInd] = max(n_multimodal_new);
        if max_multimodal > 0
            n_multimodal_new = -multimodal_inhibition .* ones(size_multimodal,1);
            n_multimodal_new(mInd) = multimodal_excitation;
        end

        %Update input buffer neurons
        if has_hipp
            if sigm
                n_buffer_new = 2.*sigmoid(.99 .* n_buffer + n_multimodal + w_buffer_hipp'*n_hipp,gain)-1; %activation
            else
                n_buffer_new = max(.99 .* n_buffer + n_multimodal +  w_buffer_hipp'*n_hipp,0); %activation
            end
        else
            if sigm
                n_buffer_new = 2.*sigmoid(.99 .* n_buffer + n_multimodal,gain)-1; %activation
            else
                n_buffer_new = max(.99 .* n_buffer + n_multimodal,0); %activation
            end
        end
        %Update pfc neurons
        if has_hipp
            if sigm
                n_pfc_new = 2.*sigmoid(w_buffer_pfc * n_buffer + 2.*w_pfc_hipp' * n_hipp,gain)-1; %activation
            else
                n_pfc_new = max(w_buffer_pfc * n_buffer + 2.*w_pfc_hipp' * n_hipp,0); %activation
            end
        else
            n_pfc_new = 2.*sigmoid(w_buffer_pfc * n_buffer,gain)-1; %activation
        end
        [max_pfc,mInd] = max(n_pfc_new);
        if max_pfc > 0
            n_pfc_new = -pfc_inhibition * ones(size_pfc,1);
            n_pfc_new(mInd) = pfc_excitation;
        end
        %Update hipp neurons
        if has_hipp
            if sigm
                n_hipp_new = 2.*sigmoid(w_pfc_hipp * n_pfc + w_buffer_hipp * n_buffer,gain)-1; %activation
            else
                n_hipp_new = max(w_pfc_hipp * n_pfc + w_buffer_hipp * n_buffer,0); %activation
            end
            %mInd = schema(i,1); % external gps signal
            [max_hipp,mInd] = max(n_hipp_new);
            if max_hipp > 0
                n_hipp_new = -hipp_inhibition * ones(size_hipp,1);
                n_hipp_new(mInd) = hipp_excitation;
            end
        end

        %Update weights
        if is_learning
            w_input_multimodal = w_input_multimodal + ojaLearningRule(lr_fast, [n_well; n_flavor], n_multimodal, w_input_multimodal, alpha_multimodal);
            w_buffer_pfc = w_buffer_pfc + ojaLearningRule(lr_slow, n_buffer, n_pfc, w_buffer_pfc, alpha_pfc);
            if has_hipp
                w_pfc_hipp = w_pfc_hipp + ojaLearningRule(lr_fast, n_pfc, n_hipp, w_pfc_hipp, alpha_hipp);
                w_buffer_hipp = w_buffer_hipp + ojaLearningRule(lr_fast, n_multimodal, n_hipp, w_buffer_hipp, alpha_hipp);
            end
        end

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

        n_multimodal = n_multimodal_new;
        n_pfc = n_pfc_new;
        n_buffer = n_buffer_new;
        if has_hipp
            n_hipp = n_hipp_new;
        end
    end

    network.w_input_multimodal = w_input_multimodal;
    network.w_buffer_pfc = w_buffer_pfc;
    network.w_pfc_hipp = w_pfc_hipp;
    network.w_buffer_hipp = w_buffer_hipp;
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