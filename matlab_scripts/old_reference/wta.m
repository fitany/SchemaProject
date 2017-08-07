% Winner take all

function wta()
    %load('Results/QR2.mat');
    load('Schema_C_07082017_2218.mat','Qr');
    QRs = Qr;
    test = QRs{243};
    if ~isempty(test)
        if ~isequal('None',test) && ~isequal([],test)
            test = test{1}
        end
        ind = str2num(test(4));
        size = 9;
        inp = .1*rand(1,size);
        if ~isempty(ind)
            inp(ind) = inp(ind)+.5;
        end
    end
    x = run_wta(size,inp);
    
    x(end)
    imagesc(x)
    
    figure;
    plot(x);
end

function sigmoid = f(x,b,a)
    % sigmoid function with threshold (b) and slope(a) parameters
    % args: x (ndarray): Input array
    %       b (float): Sigmoid threshold
    %       a (float): Sigmoid slope
    % returns: 
    %       A sigmoid computed on vector x
    sigmoid = 1 ./ (1 + exp(-((x - b) ./ a)));
end

function x = run_wta(num_neurons, inputs)
    % Main WTA routine. 
    %   Args:
    %   num_neurons (int)   : Total number of neurons
    %   inputs (ndarray)    : A numpy array containing the inputs
    %
    %   Returns:
    %   The activity (x) of the neurons (numpy array of size 
    %   sim_time x num_neurons)
    ms = 0.001;                  % ms definition :p
    dt = 10*ms;                  % Euler time step
    tf = 10;                     % total time
    tau = 100*ms;                % time-scale constant (in ms)
    sim_time = round(tf/dt);     % Total simulation time
    
    x = zeros(sim_time, num_neurons);   % Neurons initial conditions
    Winh = 1.1;                         % Inhibitory gain (global)
    Iext = inputs;                      % External inputs
    
    % Euler integration of Dynamical System (RNN)
    for t=1:sim_time
        % Excluding self-inhibition :p
        X = repmat(x(t,:), 1, num_neurons);
        X = reshape(X, num_neurons, num_neurons);
        X(logical(eye(num_neurons))) = 0;
        Iinh = sum(Winh * f(X,.5,.125),2)';
        x(t+1,:) = x(t,:) + dt/tau*(-x(t,:)-Iinh+Iext);
    end
end