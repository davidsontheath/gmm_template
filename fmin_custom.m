function best_param = fmin_custom(func, start_point)

% takes in a function func, and a vector that is
% a valid starting point for func.

% just steps around until the home point
% is the best in the neighbourhood.

size_small = 1e-5;
size_med = 1e-3;
size_big = 1e-1;

n_params = size(start_point,1);

best_param = start_point;
best_val = func(start_point);
not_done = 1;
while not_done == 1;
    
    not_done = 0;
        
    for i = 1:n_params
            step_up = best_param;
            step_up(i) = step_up(i)+size_big;
            if func(step_up) < best_val
                best_param = step_up;
                best_val = func(step_up);
                not_done = 1;
            end
            step_down = best_param;
            step_down(i) = step_down(i)-size_big;
            if func(step_down) < best_val
                best_param = step_down;
                best_val = func(step_down);
                not_done = 1;
            end
            step_up = best_param;
            step_up(i) = step_up(i)+size_med;
            if func(step_up) < best_val
                best_param = step_up;
                best_val = func(step_up);
                not_done = 1;
            end
            step_down = best_param;
            step_down(i) = step_down(i)-size_med;
            if func(step_down) < best_val
                best_param = step_down;
                best_val = func(step_down);
                not_done = 1;
            end
            step_up = best_param;
            step_up(i) = step_up(i)+size_small;
            if func(step_up) < best_val
                best_param = step_up;
                best_val = func(step_up);
                not_done = 1;
            end
            step_down = best_param;
            step_down(i) = step_down(i)-size_small;
            if func(step_down) < best_val
                best_param = step_down;
                best_val = func(step_down);
                not_done = 1;
            end
        end %i


end %while not-done


