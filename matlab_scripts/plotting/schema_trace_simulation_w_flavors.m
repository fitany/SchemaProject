%%% schema_trace_simulation_w_flavors.m
%%% Xinyun Zou
%%% Simulation of the robot's zig zag trace in different schemas with flavors in CARL lab
%%% Last updated on Aug 7, 2017
%%% Example of how to run the function:
%%% [x,y,flavor_ind,flavor_in_order]=schema_trace_simulation_w_flavors('Schema A');

function [x,y,flavor_ind,flavor_in_order]=schema_trace_simulation_w_flavors(SchemaName)
% arrays of qr positions
qr_x = [];
qr_y = [];
for i = 0:4
    for j = 0:4
        qr_x(end+1) = 60*i;
        if mod(i,2)==0
            qr_y(end+1) = 60*j;
        else
            qr_y(end+1) = 60*(4-j);
        end
    end
end

% matrix and cell to store coordinates and labels of the explored qr flavors
explored_flavor_pos = []; 
explored_flavor = {};

% two vectors to store x and y coordinates of the rest targets
waypoints_x = qr_x;
waypoints_y = qr_y;

% two vectors to store x and y coordinates of the robot trace
x = []; y = [];

% two variables to update x and y coordinates of the robot's position
updated_x = waypoints_x(1);
updated_y = waypoints_y(1);

% update robot's trace
x(end+1) = updated_x;
y(end+1) = updated_y;

% two vectors to store coordinates of the explored qr positions
prev_x = []; prev_y = [];

% two variables to update coordinates of previously explored qr position
prev_qr_x = waypoints_x(1);
prev_qr_y = waypoints_y(1);

% update the array of robot's explored qr positions
prev_x(end+1) = prev_qr_x;
prev_y(end+1) = prev_qr_y;

disp('**********Trial started**********');

% two variables to update coordinates of currently targeted qr position
dest_qr_x = waypoints_x(1);
dest_qr_y = waypoints_y(1);
waypoints_x = waypoints_x(2:end);
waypoints_y = waypoints_y(2:end);
fprintf('Targeted qr position updated to (%d, %d).\n',dest_qr_x,dest_qr_y);

ind = [];
ind(end+1) = length(x);
fprintf('ind = %d, (x,y) = (%f,%f).\n',ind(end),x(end),y(end));

flavor_ind = [];
flavor_in_order = [];

% start to plot in real-time
[explored_flavor_pos,explored_flavor,new_flavor_added]=trace_plot(qr_x,qr_y,...
    prev_x,prev_y,prev_qr_x,prev_qr_y,dest_qr_x,dest_qr_y,...
    updated_x,updated_y,explored_flavor_pos,explored_flavor,SchemaName);

if (new_flavor_added == true)
    flavor_ind(end+1) = length(x);
    flavor_in_order(end+1) = str2double(extractAfter(explored_flavor{end},2));
    fprintf('flavor_ind = %d, flavor = %s, (x,y) = (%f,%f).\n',...
        flavor_ind(end),explored_flavor{end},x(end),y(end));
end

% start simulation
isDone = false; 
isRunning = true;
while isRunning && ~isDone
    isRunning = false;
    % update targeted qr position when the current one is reached
    if(prev_qr_x == dest_qr_x && prev_qr_y == dest_qr_y)
        dest_qr_x = waypoints_x(1);
        dest_qr_y = waypoints_y(1);
        waypoints_x = waypoints_x(2:end);
        waypoints_y = waypoints_y(2:end);
        fprintf('Targeted qr position updated to (%d, %d).\n',dest_qr_x,dest_qr_y);
        if isempty(waypoints_x)
            isDone = true;
        end
    end
    % calculate the desired angle to turn 
    vectorX = dest_qr_x - prev_qr_x;
    vectorY = dest_qr_y - prev_qr_y;
    desired_angle = atan2(vectorY,vectorX);
    % use a while loop to update robot's position when it's approaching the
    % targeted position
    counter = 0;
    while (counter < ceil(4.46*norm([vectorX;vectorY])))
        updated_x = prev_qr_x + counter / 4.46 * cos(desired_angle);
        updated_y = prev_qr_y + counter / 4.46 * sin(desired_angle);
        x(end+1) = updated_x;
        y(end+1) = updated_y;
        % update the previously explored qr position if reached
        if (abs(dest_qr_x-updated_x)<0.3 && abs(dest_qr_y-updated_y)<0.3)
            ind(end+1) = length(x);
            fprintf('ind = %d, (x,y) = (%f,%f).\n',ind(end),x(end),y(end));
            prev_qr_x = dest_qr_x;
            prev_qr_y = dest_qr_y;
            prev_x(end+1) = prev_qr_x;
            prev_y(end+1) = prev_qr_y;
            isRunning = true;
        end
        % update plot in real-time
        [explored_flavor_pos,explored_flavor,new_flavor_added]=trace_plot(qr_x,qr_y,...
            prev_x,prev_y,prev_qr_x,prev_qr_y,dest_qr_x,dest_qr_y,...
            updated_x,updated_y,explored_flavor_pos,explored_flavor,SchemaName);
        if (new_flavor_added == true)
            flavor_ind(end+1) = length(x);
            flavor_in_order(end+1) = str2double(extractAfter(explored_flavor{end},2));
            fprintf('flavor_ind = %d, flavor = %s, (x,y) = (%f,%f).\n',...
                flavor_ind(end),explored_flavor{end},x(end),y(end));
        end
        % update while loop counter
        counter = counter+1;
    end
end

disp('**********Trial finished**********');

% figure(1);clf
% axis([-10 250 -10 250]);
% axis square
% hold on
% for i = 1:length(qr_x)
%     plot(qr_x(i),qr_y(i),'ro');
%     pause(0.1);
% end
% xticks(unique(qr_x));
% yticks(unique(qr_y));
% hold off

% figure(2);clf
% axis([-10 250 -10 250]);
% axis square
% hold on
% plot(qr_x,qr_y,'ro');
% for i = 1:length(x)
%     plot(x(i),y(i),'k*');
%     pause(0.001);
% end
% xticks(unique(qr_x));
% yticks(unique(qr_y));
% hold off

% for i = 1:length(x)
%     figure(3);clf
%     plot(qr_x,qr_y,'ro',x(i),y(i),'k*');
%     axis([-10 250 -10 250]);
%     axis square
%     xticks(unique(qr_x));
%     yticks(unique(qr_y));
%     pause(0.001);
% end

% figure(4);clf
% axis([-10 250 -10 250]);
% axis square
% hold on
% plot(qr_x,qr_y,'yo');
% for f = 1:length(flavorA)
%     plot(flavorA_x(f),flavorA_y(f),'o');
% end
% lgd = ['qr positions',flavorA];
% legend(lgd,'Location','bestoutside');
% xticks(unique(qr_x));
% yticks(unique(qr_y));
% hold off

end

function [explored_flavor_pos,explored_flavor,new_flavor_added]=trace_plot(qr_x,qr_y,...
    prev_x,prev_y,prev_qr_x,prev_qr_y,dest_qr_x,dest_qr_y,updated_x,updated_y,...
    explored_flavor_pos,explored_flavor,SchemaName)
% arrays of coordinates of qr flavors for Schema A
flavorA_x = [60,240,120,0,120];
flavorA_y = [240,180,0,60,120];
flavorA_pos = [flavorA_x; flavorA_y];
flavorA = {'QR01','QR02','QR03','QR04','QR05'};

% arrays of coordinates of qr flavors for Schema A swapped
flavorA_swapped_x = [60,120,0,120,240];
flavorA_swapped_y = [240,0,60,120,180];
flavorA_swapped_pos = [flavorA_swapped_x; flavorA_swapped_y];
flavorA_swapped = {'QR01','QR03','QR04','QR05','QR06'};

% arrays of coordinates of qr flavors for Schema A moved
flavorA_moved_x = [60,240,180,0,120];
flavorA_moved_y = [240,180,0,60,120];
flavorA_moved_pos = [flavorA_moved_x; flavorA_moved_y];
flavorA_moved = {'QR01','QR02','QR03','QR04','QR05'};

% arrays of coordinates of qr flavors for Schema B
flavorB_x = [60,0,180,60,240];
flavorB_y = [240,60,60,120,0];
flavorB_pos = [flavorB_x; flavorB_y];
flavorB = {'QR01','QR07','QR08','QR09','QR10'};

% arrays of coordinates of qr flavors for Schema C
flavorC_x = [60,120,240,180,0];
flavorC_y = [0,240,60,240,120];
flavorC_pos = [flavorC_x; flavorC_y];
flavorC = {'QR01','QR02','QR07','QR08','QR10'};

if (max(strcmp(SchemaName,{'Schema A','SchemaA','Schema_A'}))==true)
    flavor_pos = flavorA_pos;
    flavor = flavorA;
elseif (max(strcmp(SchemaName,{'Schema B','SchemaB','Schema_B'}))==true)
    flavor_pos = flavorB_pos;
    flavor = flavorB;
elseif (max(strcmp(SchemaName,{'Schema C','SchemaC','Schema_C'}))==true)
    flavor_pos = flavorC_pos;
    flavor = flavorC;
elseif (max(strcmp(SchemaName,{'Schema A swapped','SchemaA swapped','Schema_A_swapped'}))==true)
    flavor_pos = flavorA_swapped_pos;
    flavor = flavorA_swapped;
elseif (max(strcmp(SchemaName,{'Schema A moved','SchemaA moved','Schema_A_moved'}))==true)
    flavor_pos = flavorA_moved_pos;
    flavor = flavorA_moved;
else
    error('Please enter a valid SchemaName.');
end

new_flavor_added = false;
[Lia,Locb]=ismember([prev_qr_x,prev_qr_y],flavor_pos','rows');
if (Lia == true)
    if ~isempty(explored_flavor_pos)
        [exp_Lia,exp_Locb]=ismember([prev_qr_x,prev_qr_y],explored_flavor_pos','rows');
        if (exp_Lia == false)
            explored_flavor_pos(:,end+1)=[prev_qr_x;prev_qr_y];
            explored_flavor{end+1}=flavor{Locb};
            new_flavor_added = true;
        end
    else
        explored_flavor_pos(:,end+1)=[prev_qr_x;prev_qr_y];
        explored_flavor{end+1}=flavor{Locb};
        new_flavor_added = true;
    end
end
% figure(5);clf
% axis([-10 250 -10 250]);
% axis square;
% hold on
% plot(qr_x,qr_y,'ro',prev_x,prev_y,'go',...
%     dest_qr_x,dest_qr_y,'bo',updated_x,updated_y,'k*');
% lgd = {'qr positions','explored places','targeted place','robot position'};
% if ~isempty(explored_flavor)
%     for f = 1:length(explored_flavor)
%         plot(explored_flavor_pos(1,f),explored_flavor_pos(2,f),'o');
%     end
%     lgd = [lgd,explored_flavor];
% end
% legend(lgd,'Location','bestoutside');
% xticks(unique(qr_x));
% yticks(unique(qr_y));
% title(SchemaName);
% hold off
% pause(1e-8);
end
