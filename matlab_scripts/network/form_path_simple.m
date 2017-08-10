%%% form_path_simple.m
%%% Author: Xinyun Zou
%%% Find out the optimal path involving upto 25 available qr positions from a
%%% starting position to all flavors detected earlier in the schema.
%%% Last updated on Aug 10, 2017
%%% Example of how to run the function:
%%% total_path = form_path_simple(startPt,simple_qr_pos,pos,flavor_ind,flavor_in_order);

function total_path = form_path_simple(startPt,simple_qr_pos,pos,flavor_ind,flavor_in_order)
[flavor_reorder,flavor_qr_pos_reorder,flavor_pos_reorder] = test_distance(startPt,pos,flavor_ind,flavor_in_order);
% qr_x = [];
% qr_y = [];
% for i = 0:4
%     for j = 0:4
%         qr_x(end+1) = 60*i;
%         if mod(i,2)==0
%             qr_y(end+1) = 60*j;
%         else
%             qr_y(end+1) = 60*(4-j);
%         end
%     end
% end
%qr_pos = [qr_x;qr_y];
flavors_pos_to_go = flavor_qr_pos_reorder;
disp(flavor_qr_pos_reorder);
endPt = flavor_qr_pos_reorder(:,end);
if isequal(startPt,flavor_qr_pos_reorder(:,1))
    total_path = [];
else
    total_path = startPt;
end
while ~isempty(flavors_pos_to_go)
    toGoPt = flavors_pos_to_go(:,1);
    pts_to_pass = form_path_between_two_pts(startPt,toGoPt,simple_qr_pos);
    total_path = [total_path, pts_to_pass];
    flavors_pos_to_go(:,1) = [];
    startPt = toGoPt;
end

end

function [flavor_reorder,flavor_qr_pos_reorder,flavor_pos_reorder] = test_distance(startPt,pos,flavor_ind,flavor_in_order)
if ~isequal(size(startPt),[2,1])
    error('Please enter a startPt input of size (2,1)');
end
%[x,y,flavor_ind,flavor_in_order]=schema_trace_simulation_w_flavors('Schema A');
flavor_pos = pos(:,flavor_ind); %flavor_pos = [x(flavor_ind);y(flavor_ind)];
flavor_qr_pos = round(flavor_pos);
flavor_waypoints = flavor_qr_pos;
flavor_qr_pos_reorder = [];
flavor_pos_reorder = [];
flavor_reorder = [];
[Lia, Locb] = ismember(startPt',flavor_qr_pos','rows');
if Lia == true
    flavor_qr_pos_reorder(:,end+1) = startPt;
    flavor_reorder(end+1) = flavor_in_order(Locb);
    flavor_pos_reorder(:,end+1) = pos(:,flavor_ind(Locb));
    %flavor_pos_reorder(:,end+1) = [x(flavor_ind(Locb));y(flavor_ind(Locb))];
end
while ~isempty(flavor_waypoints)
    [startPt,flavor_waypoints] = findNearestPoint(flavor_waypoints,startPt);
    flavor_qr_pos_reorder(:,end+1)=startPt;
    [Lia, Locb] = ismember(startPt',flavor_qr_pos','rows');
    flavor_reorder(end+1) = flavor_in_order(Locb);
    flavor_pos_reorder(:,end+1) = pos(:,flavor_ind(Locb));
    %flavor_pos_reorder(:,end+1) = [x(flavor_ind(Locb));y(flavor_ind(Locb))];
end
end

function pts_to_pass = form_path_between_two_pts(startPt,endPt,qr_pos)
fprintf('startPt is (%d,%d), and endPt is (%d,%d).\n',...
    startPt(1),startPt(2),endPt(1),endPt(2));
pts_to_pass = [];
backup_pts = [];
if norm(startPt-endPt)>60*sqrt(2)
    for i = 1:size(qr_pos,2)
        d = point_to_line(qr_pos(:,i), startPt, endPt);
        if (d < 30*sqrt(2) && d >= 0 && ~isequal(qr_pos(:,i),startPt) && ~isequal(qr_pos(:,i),endPt))
            angle1 = abs(angle_between_lines(qr_pos(:,i)-endPt,startPt-endPt));
            angle2 = abs(angle_between_lines(qr_pos(:,i)-startPt,endPt-startPt));
            if (angle1 < pi/2 && angle2 < pi/2)
                backup_pts(:,end+1)=qr_pos(:,i);
            end
        end
    end
end

while ~isempty(backup_pts)
    [nextPt,backup_pts] = findNearestPoint(backup_pts,startPt);
    pts_to_pass(:,end+1) = nextPt;
end

pts_to_pass(:,end+1) = endPt;
pts_to_pass = shrink_path(pts_to_pass);
disp(pts_to_pass);

end

function d = point_to_line(pt, v1, v2)
a = v1 - v2;
b = pt - v2;
d = sqrt(norm(b)^2 - dot(b,a)^2/norm(a)^2);
end

function angle = angle_between_lines(v1,v2)
v1 = [v1;0];
v2 = [v2;0];
angle = atan2(norm(cross(v1,v2)),dot(v1,v2));
end

function [nextPt,flavor_waypoints] = findNearestPoint(flavor_waypoints,prevPt)
diff = flavor_waypoints - prevPt;
distance = zeros(1,size(flavor_waypoints,2));
for i = 1:size(flavor_waypoints,2)
    distance(i) = norm(diff(:,i));
end
[m,ind]=min(distance);
nextPt = flavor_waypoints(:,ind);
flavor_waypoints(:,ind) = [];
end

function pts_to_pass = shrink_path(pts_to_pass)
while size(pts_to_pass,2) >= 3
%     disp('#######');
%     disp(pts_to_pass);
%     disp('#######');
    orig_length = size(pts_to_pass,2);
    for i = size(pts_to_pass,2):-1:3
        dist1 = norm(pts_to_pass(:,i)-pts_to_pass(:,i-1));
        dist2 = norm(pts_to_pass(:,i)-pts_to_pass(:,i-2));
        if (dist1<=60*sqrt(2) && dist2<=60*sqrt(2))
            pts_to_pass(:,i-1)=[];
            continue
        end
    end
    if (size(pts_to_pass,2) == orig_length)
        break
    end
end
end