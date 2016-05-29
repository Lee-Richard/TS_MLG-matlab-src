%function TSMG
%
%inputs:
%simulation time: simulation time in seconds
%lambda: big flow arrival rate
%miu: inverse of flow lifetime
%
%outputs:
%
%

function [data_of_interest_m] = tsmg(simulation_time,lambda,miu,max_storage_times,temporal_cost,node_storage,wavelength)



%function TSMG
%
%inputs:
%simulation time: simulation time in seconds
%lambda: big flow arrival rate
%miu: inverse of flow lifetime
%
%outputs:
%
%

%simulation parameters
COST = temporal_cost;    %存储代价
SIMULATIONTIME = simulation_time;    %仿真时间
SLOTPERSECOND=50;    %时隙
TIMES = max_storage_times;    %最大存储次数，映射为多层图层数
STORAGE = node_storage*SLOTPERSECOND;    %每个节点的存储大小
WAVELENGTH = wavelength;    %每条链路上的波长数

%constants
SHUTDOWN_PERIOD = 0;    %仿真关闭时间，为0表示不关闭
%Unify bandwidth values with flow rate
flow_rate = 1;    %流速率，乘以持续时间即为流大小
applied_bandwidth = 1;    %单个业务请求占用的带宽

%Run time variables
flow_id = 0;     %流id，递增，唯一

%data of interest
total_delay = 0;  %总延迟
bandwidth_utilization=0;  %链路利用率
number_of_rejected_flow = 0;    %拒绝的流
number_of_generated_flow = 0;    %总的流
number_of_blocking_flow = 0;    %缓存的流
capacity_of_rejected_flow = 0;    %拒绝的流大小
capacity_of_all_flow = 0;     %流的总大小
used_bandwidth = 0;    %已经被用过的总带宽
hop_count = 0;    %？
max_storage_utilization = zeros(1,14);    %存储最大利用率
number_of_storaged_flow=0;    %被存储的流的数量
max_layers = 0;

%%auxiliary tables and parameters
node_numbers=14;
flow_table = [];    %记录活动流的信息
zero_auxiliary_vector = zeros(1,node_numbers);
inf_auxiliary_vector = ones(1,node_numbers)*inf;
layer_interval_table = 0;    %行向量，记录第i层与第一层之间的时间间隔
storage_table = diag(ones(1,node_numbers)*STORAGE); 
storage_cost_table=diag(ones(1,node_numbers))./COST;    %记录节点存储使用的信息

adjacency_matrix = [inf,1,inf,inf,inf,inf,1,inf,inf,inf,inf,inf,inf,1;
                    1,inf,inf,inf,1,inf,inf,inf,inf,inf,inf,inf,inf,1;
                    inf,inf,inf,1,inf,inf,inf,inf,1,inf,inf,inf,inf,1;
                    inf,inf,1,inf,1,1,inf,inf,inf,inf,inf,inf,inf,inf;
                    inf,1,inf,1,inf,inf,inf,1,inf,inf,inf,inf,inf,inf;
                    inf,inf,inf,1,inf,inf,1,inf,inf,inf,inf,inf,inf,inf;
                    1,inf,inf,inf,inf,1,inf,inf,inf,1,inf,inf,inf,inf;
                    inf,inf,inf,inf,1,inf,inf,inf,inf,1,inf,inf,inf,inf;
                    inf,inf,1,inf,inf,inf,inf,inf,inf,inf,1,inf,inf,inf;
                    inf,inf,inf,inf,inf,inf,1,1,inf,inf,1,inf,1,inf;
                    inf,inf,inf,inf,inf,inf,inf,inf,1,1,inf,1,inf,inf;
                    inf,inf,inf,inf,1,inf,inf,inf,inf,inf,1,inf,1,inf;
                    inf,inf,inf,inf,inf,inf,inf,inf,1,1,inf,1,inf,inf;
                    1,1,1,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf,inf
                    ];                         %  14点的NSFNET的邻接矩阵
% adjacency_matrix = [inf,1,1;
%                     1,inf,1;
%                     1,1,inf];
adjacency_auxiliary_matrix = 1./storage_cost_table;

unit_bandwidth_table =   [0,1,0,0,0,0,1,0,0,0,0,0,0,1;
                          1,0,0,0,1,0,0,0,0,0,0,0,0,1;
                          0,0,0,1,0,0,0,0,1,0,0,0,0,1;
                          0,0,1,0,1,1,0,0,0,0,0,0,0,0;
                          0,1,0,1,0,0,0,1,0,0,0,0,0,0;
                          0,0,0,1,0,0,1,0,0,0,0,0,0,0;
                          1,0,0,0,0,1,0,0,0,1,0,0,0,0;
                          0,0,0,0,1,0,0,0,0,1,0,0,0,0;
                          0,0,1,0,0,0,0,0,0,0,1,0,0,0;
                          0,0,0,0,0,0,1,1,0,0,1,0,1,0;
                          0,0,0,0,0,0,0,0,1,1,0,1,0,0;
                          0,0,0,0,1,0,0,0,0,0,1,0,1,0;
                          0,0,0,0,0,0,0,0,1,1,0,1,0,0;
                          1,1,1,0,0,0,0,0,0,0,0,0,0,0
                          ];                            %单位带宽的矩阵
% unit_bandwidth_table = [0,1,1;
%                         1,0,1;
%                         1,1,0];
bandwidth_per_link = WAVELENGTH;                           %每条空间链路上的波长
bandwidth_table = unit_bandwidth_table * bandwidth_per_link;       %单层图的链路可用带宽表
tsml_bandwidth_table = [bandwidth_table, storage_table];      %空间域和时间域共同组成了多层图
tsml_table = [adjacency_matrix,adjacency_auxiliary_matrix];  
total_bandwidth=sum(sum(bandwidth_table));    %  网络总带宽

fprintf('the init tsml table\n');
disp(tsml_table);
fprintf('the init tsml_bandwidth_table\n');
disp(tsml_bandwidth_table);

%main loop
% fid=fopen('tsml.log','wt');
% fprintf(fid,'\ntsml table at 0(beginning):\n');
% fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',tsml_table');
% fprintf(fid,'\ntsml bandwidth table at 0(beginning):\n');
% fprintf(fid,'%d %d %d %d %d %d %d %d %d\n',tsml_bandwidth_table');
% fprintf(fid,'\nlayer table at 0(beginning):\n');
% fprintf(fid,'%d %d %d %d %d\n',layer_table');
% fprintf(fid,'------------------------------------------------------------------------------------------------\n');
% disp(tsml_table);
% disp(tsml_bandwidth_table);


time_to_wait = uint32(exprnd(1/lambda)*SLOTPERSECOND);
%fprintf(fid,'\ntime to wait = %d\n',time_to_wait);
% rand('state',0);
for ite = 1:SIMULATIONTIME*SLOTPERSECOND
    %fprintf(fid,'Round %d start:\n',ite);
%     time_to_wait = uint32(50 * exprnd(1/lambda));
%     time_to_wait = time_to_wait-1;
fprintf('ite=%d\n',ite);
    if time_to_wait == 0                     % 生成新流
        time_to_wait = uint32(exprnd(1/lambda)*SLOTPERSECOND);
        %fprintf(fid,'\ntime to wait = %d\n',time_to_wait);
        if ite< (SIMULATIONTIME - SHUTDOWN_PERIOD)*SLOTPERSECOND
            duration = ceil(exprnd(1/miu)*SLOTPERSECOND);
            source = randi(node_numbers);
            destination = randi(node_numbers);
            while(destination == source)
                destination = randi(node_numbers);
            end
            flow_id = flow_id + 1;
            file_size = duration * flow_rate;      %文件大小
            capacity_of_all_flow = capacity_of_all_flow + file_size;
            %fprintf(fid,'New flow generated at round %d: id=%d, source=%d, destination=%d,duration=%d\n',ite,flow_id, source, destination,duration);
            new_flow = [flow_id, duration, file_size, source, destination];
            flow_table=[flow_table;new_flow];
            number_of_generated_flow = number_of_generated_flow + 1;
            fprintf('new flow: source=%d, destination=%d, duration=%d\n',source,destination,duration);
        end
    
        %构造辅助图 
        auxiliary_tsml_table = tsml_table;
        layers = size(tsml_table,1)/node_numbers;    %计算当前时移多层图的层数
        
        if layers> max_layers
            max_layers = layers;
        end
        disp(layers);
        disp(max_layers);
%         
%         %去掉剩余带宽小于请求带宽(默认为1)的空间链路
%         if layers ~= 1
%             %查找tsml_bandwidth_table中，所有带宽小于flow_rate的位置被移除
%         end
%         
%         %去掉存储容量小于业务请求文件大小的时域链路
        if layers ~= 1
            %查找tsml_bandwidth_table中，所有存储小于fize_size的位置被移除
            for i = 1:size(tsml_bandwidth_table,1)
                if tsml_bandwidth_table(i,i+node_numbers) < file_size
                    auxiliary_tsml_table(i,i+node_numbers) = inf;
                end
            end
        end
        
        %去掉空闲时间小于业务请求持续时间的时间链路
        if layers ~= 1
            %查找空闲时间小于duration的链路，移除
            [row, col] = find(tsml_bandwidth_table == 1);
            idle_threshold = [row, col];
            for i = 1:size(idle_threshold,1)
                idle_row = idle_threshold(i,1);
                idle_col = idle_threshold(i,2);
                layer_idle = floor((idle_row-1)/node_numbers) + 1;
                while (idle_row <= size(tsml_bandwidth_table,1) && idle_col <= size(tsml_bandwidth_table,2))    %查找link的首次忙碌时间
                    idle_row = idle_row + node_numbers;
                    idle_col = idle_col + node_numbers;
                    if idle_row > size(tsml_bandwidth_table,1) 
                        break;
                    end
                    if tsml_bandwidth_table(idle_row,idle_col) == 0
                        %记录下空闲时间
                        layer_busy = floor((idle_row-1)/node_numbers) + 1;
                        idle_time = layer_interval_table(1,layer_busy) - layer_interval_table(1,layer_idle);
                        %更新辅助图
                        if idle_time < duration
                            auxiliary_tsml_table(idle_threshold(i,1),idle_threshold(i,2)) = inf;
                        end
                        break;
                    end
                end
            end
        end

        %在刚才构造的辅助图中查找路径
        i=1;

        while(i <= layers)
            [distance,path] = dijkstra(auxiliary_tsml_table,source,destination+(i-1)*node_numbers);
            if distance == inf
                i = i+1;
                continue;
            else
                break;
            end
        end

        fprintf('path = ');
        disp(path);
        
        %if path exists, update tsmlg
        if distance < inf
            %在tsml_bandwidth_table中减去被占用的带宽，包含存储
%             for i=1:size(path,2)-1
%                 tsml_bandwidth_table(path(i),path(i+1)) = tsml_bandwidth_table(path(i),path(i+1)) - applied_bandwidth;
%                 if tsml_bandwidth_table(path(i),path(i+1)) == 0     %同步更新tsml_table
%                     tsml_table(path(i), path(i+1)) = inf;
%                 end
%             end
            
            %%%%%%更新时移多层图，用到了几层图就添加几层图,遍历path,每一层加一个图。新流的持续时间，和各个层的间隔比较
            %如果加层后，总层数超过了TIMES,就直接break
            storage_switch = 1;    %初始化：存→路or存→存，为了加层
            % index1 = 1;     
            index_end = 1;    %每次开始前重置是必要的
            storage_index_start = -1;
            storage_index_end = 0;
            layer_overlay_flag = 0;
            %j=1;    %保证interval_table和tsml_table的层数一致
            temp_layer_interval_table = layer_interval_table;
            reserve_layer_interval_table = layer_interval_table;
            reserve_tsml_table = tsml_table;
            reserve_tsml_bandwidth_table = tsml_bandwidth_table;
            j=1;
            while (j < size(path,2)-1)    %每一轮循环过两个点，更新一轮多层图和interval_table(可选)
                if mod(path(j),node_numbers) == mod(path(j+1),node_numbers)    %时域路径，存
                    layer = floor((path(j)-1)/node_numbers) + 1;
                    if storage_switch == 0    %开始总结存储 路→存
                        storage_switch = 1;
                        layer = floor((path(j)-1)/node_numbers) + 1;
                        %记录下index1
                        %storage_index_start = floor(path(j)/(node_numbers+1)) + 1;    %有问题吗?
                        storage_index_start = find (temp_layer_interval_table == layer_interval_table(layer));
                        j = j+1;
                        continue;
                    else    %继续总结存储   存→存
                        if j == 1
                            storage_index_start = find (temp_layer_interval_table == layer_interval_table(layer));
                            index_end = floor((path(j+1)-1)/node_numbers)+1;
                            disp('aaaaaa');
                            storage_switch = 1;
                        end
                        j = j+1;
                        continue;
                    end
                else    %空间路径，带宽         
                    if storage_switch ==1    %存→路，可能加层，记录index2，先更新index1和index2，更新存储，
                        storage_switch = 0;
                        %更新index1和index2；
                        layer = floor((path(j)-1)/node_numbers) + 1;
                        interval_start = max(temp_layer_interval_table(index_end),layer_interval_table(layer));   %取上一层的index_start+duration or 上层的index_end之间的较大值？
                        interval_end = interval_start + duration;    
                        if any(layer_interval_table == interval_end)    %有重叠
                            %do sth,不加层
                            index_start = find(temp_layer_interval_table == interval_start);    %
                            index_end = find(temp_layer_interval_table == interval_end);
                            storage_index_end = index_end;    %新加的
                             %减带宽资源
                            for k=index_start:index_end-1   %层
                                tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) = ...
                                    tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) - applied_bandwidth;
                                if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) == 0     %同步更新tsml_table
                                    tsml_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) = inf;
                                end
                                if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) < 0     %异常处理
                                    disp('268');
                                    return;
                                end
                            end

                            %减存储资源
                            if storage_index_start ~= -1    %排除初始化的情形
                               for k=storage_index_start:storage_index_end-1   %层
                                    tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) =...
                                        tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) - file_size;
                                    if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) == 0     %同步更新tsml_table
                                        tsml_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) = inf;
                                    end
                                    if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) < 0     %异常处理
                                        disp('282');
                                        return;
                                    end
                               end 
                            end          
                        else    %无重叠层
            %                temp_layer_interval_table = [temp_layer_interval_table, interval_end];
                            %判断是否超出了层
                            if size(tsml_table,1) == TIMES*node_numbers;
                                tsml_bandwidth_table = reserve_tsml_bandwidth_table;
                                tsml_table = reserve_tsml_table;
                                temp_layer_interval_table = reserve_layer_interval_table;
                                disp('reject the flow');
                                number_of_rejected_flow = number_of_rejected_flow + 1;
                                break;
                            end
                            temp_layer_interval_table = [temp_layer_interval_table, interval_end];
                            temp_layer_interval_table = sort(temp_layer_interval_table);
                            index_start = find(temp_layer_interval_table == interval_start);    %
                            index_end = find(temp_layer_interval_table == interval_end);
                            storage_index_end = index_end;
                            %复制index_end-1层到index_end层，然后，index之上的所有行，尾部加零;index（包含）之下的所有行，前面加零;
                            temp_tsml_bandwidth_table = [];
                            temp_tsml_table = [];
                            %disp(index_end);
                            for k=(index_end-1)*node_numbers+1:index_end*node_numbers    %???? index=3, k=7
                                tsml_bandwidth_table = [tsml_bandwidth_table(1:k-1,:);tsml_bandwidth_table(k-node_numbers,:);tsml_bandwidth_table(k:end,:)];
                                tsml_table = [tsml_table(1:k-1,:);tsml_table(k-node_numbers,:);tsml_table(k:end,:)];
                            end
                            %每行加零,后头加
                            for k=1:(index_end-1)*node_numbers
                                temp_tsml_bandwidth_table = [temp_tsml_bandwidth_table,tsml_bandwidth_table(k,:),zero_auxiliary_vector];     %语法有错误 reshape函数？
                                temp_tsml_table = [temp_tsml_table,tsml_table(k,:),inf_auxiliary_vector]; %语法有错误
                            end
                            %每行加零，前头加
                            for k=(index_end-1)*node_numbers+1:size(tsml_bandwidth_table,1)
                                temp_tsml_bandwidth_table = [temp_tsml_bandwidth_table,zero_auxiliary_vector,tsml_bandwidth_table(k,:)]; %语法有错误
                                temp_tsml_table = [temp_tsml_table,inf_auxiliary_vector,tsml_table(k,:)]; %语法有错误
                            end
                            tsml_bandwidth_table = reshape(temp_tsml_bandwidth_table, (size(temp_layer_interval_table,2)+1)*node_numbers,size(temp_layer_interval_table,2)*node_numbers);
                            tsml_bandwidth_table = tsml_bandwidth_table.';
                            tsml_table = reshape(temp_tsml_table, (size(temp_layer_interval_table,2)+1)*node_numbers,size(temp_layer_interval_table,2)*node_numbers);
                            tsml_table = tsml_table.';

                            %disp(tsml_bandwidth_table);

                            %减带宽资源
                          for k=index_start:index_end-1   %层
                                tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) = ...
                                    tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) - applied_bandwidth;
                              if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) < 0  %层重叠
                                  auxiliary_tsml_table(path(j),path(j+1)) = inf;
                                  layer_overlay_flag = 1;
                                  break;
                              end
                              if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) == 0     %同步更新tsml_table
                                    tsml_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) = inf;
                              end
%                                 if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) < 0     %异常处理
%                                     disp('335');
%                                     return;
%                                 end
                          end
                          if layer_overlay_flag == 1    %重叠
                              %找路径，初始化，continue
                              layer_overlay_flag = 0;
                                i=1;
                                while(i <= layers)
                                    [distance,path] = dijkstra(auxiliary_tsml_table,source,destination+(i-1)*node_numbers);
                                    if distance == inf
                                        i = i+1;
                                        continue;
                                    else
                                        break;
                                    end
                                end
                                if distance < inf
                                    j=1;
                                    storage_switch = 1;    %初始化：存→路or存→存，为了加层
                                    % index1 = 1;     
                                    index_end = 1;    %每次开始前重置是必要的
                                    storage_index_start = -1;
                                    storage_index_end = 0;
                                    layer_overlay_flag = 0;
                                    %j=1;    %保证interval_table和tsml_table的层数一致
                                    temp_layer_interval_table = layer_interval_table;
                                    layer_interval_table=reserve_layer_interval_table ;
                                    tsml_table =reserve_tsml_table;
                                    tsml_bandwidth_table = reserve_tsml_bandwidth_table;
                                    continue;
                                else
                                    disp('reject the flow');
                                    number_of_rejected_flow = number_of_rejected_flow +1;
                                    break;
                                end
                          end
                           
                           %减存储资源
                           if storage_index_start ~= -1    %排除初始化的情形,初始化默认上层为存，不减存储的存
                               for k=storage_index_start:storage_index_end-1   %层
                                    tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) =...
                                        tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) - file_size;
                                    if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) == 0     %同步更新tsml_table
                                        tsml_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) = inf;
                                    end
                                    if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j)-1,node_numbers)+1+(k)*node_numbers) < 0     %异常处理
                                        disp('349');
                                        return;
                                    end
                               end 
                           end 
                        end  
                    else    %路→路，正常,续用之前的index1和index2，不加层，减资源  
                        for k=index_start:index_end-1   %层
                            tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) = ...
                                tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) - applied_bandwidth;
                            if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) < 0
                                auxiliary_tsml_table(path(j),path(j+1)) = inf;
                                layer_overlay_flag = 1;
                                break;
                            end
                            if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) == 0     %同步更新tsml_table
                                tsml_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) = inf;
                            end
%                             if tsml_bandwidth_table(mod(path(j)-1,node_numbers)+1+(k-1)*node_numbers,mod(path(j+1)-1,node_numbers)+1+(k-1)*node_numbers) < 0     %异常处理
%                                 disp('363');
%                                 return;
%                             end
                        end
                        if layer_overlay_flag == 1    %重叠
                          %找路径，初始化，continue
                          layer_overlay_flag = 0;
                            i=1;
                            while(i <= layers)
                                [distance,path] = dijkstra(auxiliary_tsml_table,source,destination+(i-1)*node_numbers);
                                if distance == inf
                                    i = i+1;
                                    continue;
                                else
                                    break;
                                end
                            end
                            if distance < inf
                                j=1;
                                storage_switch = 1;    %初始化：存→路or存→存，为了加层
                                % index1 = 1;     
                                index_end = 1;    %每次开始前重置是必要的
                                storage_index_start = -1;
                                storage_index_end = 0;
                                layer_overlay_flag = 0;
                                %j=1;    %保证interval_table和tsml_table的层数一致
                                temp_layer_interval_table = layer_interval_table;
                                layer_interval_table=reserve_layer_interval_table ;
                                tsml_table =reserve_tsml_table;
                                tsml_bandwidth_table = reserve_tsml_bandwidth_table;
                                continue;
                            else
                                disp('reject the flow');
                                number_of_rejected_flow = number_of_rejected_flow + 1;
                                break;
                            end
                        end
                    end

                end
                %disp(tsml_bandwidth_table);   
%                 fprintf('j=%d\n',j);
%                 disp(tsml_bandwidth_table);
                j=j+1;
            end
            layer_interval_table = temp_layer_interval_table;
%             fprintf('before transmittion, the layer interval table is\n');
%             disp(layer_interval_table);
            %处理存储资源的占用
            %total_delay = total_delay+delay;
            %新流用在时域上分成几段传输完成的
            %segments = path(1,end)/node_numbers+1;   %？    
        %path not exists
        else
            %reject the flow
            disp('reject the flow');
            number_of_rejected_flow = number_of_rejected_flow +1;
        end
    else
        %fprintf(fid,'No new flow generated at round %d due to flow interval\n',ite);
        
    end
    time_to_wait = time_to_wait -1;
%     disp(time_to_wait);
%     disp(ite);
%判断是否有流即将传送完毕,所有层的interval减1，是否可行？减一层，并更新layer_interval_table？
    if size(layer_interval_table,2) > 1
        for i=2:size(layer_interval_table,2)
            layer_interval_table(1,i) = layer_interval_table(1,i) - 1;
        end
        if layer_interval_table(1,2) == 0  %如果有传输完毕的情况，删除一层(第二层)
            layer_interval_table(:,1) = [];
            tsml_bandwidth_table(1:node_numbers,:) = [];
            tsml_bandwidth_table(:,1:node_numbers) = [];
            tsml_table(1:node_numbers,:) = [];
            tsml_table(:,1:node_numbers) = [];
        end
    end
%统计结果

% fprintf('at the end of %d ite, the tsml_table\n', ite);
% disp(tsml_table);
% fprintf('at the end of %d ite, the tsml_bandwidth_table\n', ite);
% disp(tsml_bandwidth_table);
% fprintf('at the end of %d ite, the layer_interval_table\n', ite);
% disp(layer_interval_table);
% fprintf('\n\n\n');
    testb = find(tsml_bandwidth_table<0);
    if size(testb,1) ~=0 
        return;
    end
end

blocking_ratio = number_of_rejected_flow/number_of_generated_flow;



data_of_interest_m = [blocking_ratio];
