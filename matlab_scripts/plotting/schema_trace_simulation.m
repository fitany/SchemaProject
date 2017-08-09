%%% schema_trace_simulation.m
%%% Xinyun Zou
%%% Simulation of the robot's zig zag trace in schema without flavor in CARL lab
%%% Last updated on Aug 3, 2017

clear;clc;

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
fprintf('ind = %d, (x,y) = (%f,%f).\n',length(x),x(end),y(end));

% % start to plot in real-time
% figure(4);clf
% plot(qr_x,qr_y,'ro',prev_x,prev_y,'go',...
%     dest_qr_x,dest_qr_y,'bo',updated_x,updated_y,'k*');
% axis([-10 250 -10 250]);
% axis square;
% lgd = {'qr positions','explored places','targeted place','robot position'};
% legend(lgd,'Location','bestoutside');
% title('Position Only (with no flavor)');
% pause(1e-5);

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
            disp('**********Trial finished**********');
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
    while (counter < ceil(4.46*60))
        updated_x = prev_qr_x + counter / 4.46 * cos(desired_angle);
        updated_y = prev_qr_y + counter / 4.46 * sin(desired_angle);
        x(end+1) = updated_x;
        y(end+1) = updated_y;
        % update the previously explored qr position if reached
        if (abs(dest_qr_x-updated_x)<0.3 && abs(dest_qr_y-updated_y)<0.3)
            ind(end+1) = length(x);
            fprintf('ind = %d, (x,y) = (%f,%f).\n',length(x),x(end),y(end));
            prev_qr_x = dest_qr_x;
            prev_qr_y = dest_qr_y;
            prev_x(end+1) = prev_qr_x;
            prev_y(end+1) = prev_qr_y;
            isRunning = true;
        end
%         % update plot in real-time
%         figure(4);clf
%         plot(qr_x,qr_y,'ro',prev_x,prev_y,'go',...
%             dest_qr_x,dest_qr_y,'bo',updated_x,updated_y,'k*');
%         axis([-10 250 -10 250]);
%         axis square;
%         lgd = {'qr positions','explored places','targeted place','robot position'};
%         legend(lgd,'Location','bestoutside');
%         xticks(unique(qr_x));
%         yticks(unique(qr_y));
%         title('Position Only (with no flavor)');
%         pause(1e-5);
        % update while loop counter
        counter = counter+1;
    end
end

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

figure(2);clf
axis([-10 250 -10 250]);
axis square
hold on
plot(qr_x,qr_y,'ro');
for i = 1:length(x)
    plot(x(i),y(i),'k*');
    pause(0.001);
end
xticks(unique(qr_x));
yticks(unique(qr_y));
hold off

% for i = 1:length(x)
%     figure(3);clf
%     plot(qr_x,qr_y,'ro',x(i),y(i),'k*');
%     axis([-10 250 -10 250]);
%     axis square
%     xticks(unique(qr_x));
%     yticks(unique(qr_y));
%     pause(0.001);
% end
