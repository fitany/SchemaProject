%%% schema_data_plot.m
%%% Xinyun Zou
%%% Plotting of the real trace data from the robot's odometry.
%%% Last updated on Aug 3, 2017

clear;clc;

if ispc
    filename = '..\..\training_data\2017-08-03 132845.28_SchemaA_moved.csv';
else
    filename = '../../training_data/2017-08-03 132845.28_SchemaA_moved.csv';
end
M = csvread(filename,1,0);
time = M(:,1)-M(1,1);
x = M(:,2); % current x position, might be wrongly calculated for earlier data files
y = M(:,3); % current x position, might be wrongly calculated for earlier data files
prev_x = M(:,4);
prev_y = M(:,5);
dest_x = M(:,6);
dest_y = M(:,7);
distance = M(:,8);
curr_angle = M(:,9);
desired_angle = M(:,10);
angle_to_turn = M(:,11);
flavor = M(:,12); % 0 if no qr flavor detected at that time
updated_x = prev_x + distance .* cos(desired_angle/180*pi); % to replace wrongly calculated x for earlier data files
updated_y = prev_y + distance .* sin(desired_angle/180*pi); % to replace wrongly calculated y for earlier data files

% Note: for old data files with only 7 columns and no column names on the first row:
% M = csvread(filename);
% time = M(:,1)-M(1,1);
% x =  M(:,2); % calculated wrongly
% y = M(:,3); % calculated wrongly
% qr_x = M(:,4);
% qr_y = M(:,5);
% distance = M(:,6);
% angle = M(:,7);


% figure(1);
% %scsz = get(0,'ScreenSize');
% %set(gcf,'Position',scsz);
% %plot(x,y,'b-',x,y,'k*');
% plot(updated_x,updated_y,'b-',updated_x,updated_y,'k*');
% hold on 
% plot(prev_x,prev_y,'ro');
% hold off
% axis([-10 250 -10 250]);
% axis square;

figure(2);
axis([-10 250 -10 250]);
axis square;
hold on
plot(prev_x,prev_y,'ro');
for i = 1:length(updated_x)
    plot(prev_x(i),prev_y(i),'go',updated_x(i),updated_y(i),'k*');
    pause(0.00001);
end
hold off
