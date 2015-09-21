clear all;

load('./may_first_week.mat');

% Calculations

arr_delay = [];
arr_delay_count = 0;
for i = 1:size(arr_data)
    if (arr_data(i) > 0)
        arr_delay_count = arr_delay_count+1;
        arr_delay(arr_delay_count) = arr_data(i);
    end
end

dep_delay = [];
dep_delay_count = 0;
for i = 1:size(dep_data)
    if (dep_data(i) > 0)
        dep_delay_count = dep_delay_count+1;
        dep_delay(dep_delay_count) = dep_data(i);
    end
end

week_delay = [];
week_delay_count = 0;
for i = 1:size(week_data)
    if (week_data(i) > 0)
        week_delay_count = week_delay_count+1;
        week_delay(week_delay_count) = week_data(i);
    end
end

% Plots

