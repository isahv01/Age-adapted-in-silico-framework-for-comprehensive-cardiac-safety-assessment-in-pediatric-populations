function T_end = visualize_TRA_method(ecg_data, fs, T_peak, R_peak, visual)
    % VISUALIZE_TRA_METHOD
    % --------------------------------------------------------------------
    % INPUTS:
    %   ecg_signal - data from one lead
    %   fs         - sample frequency (Hz)
    %   T_peak     - peak from T wave (sample)
    %   R_peak     - peak from R wave (sample)
    % --------------------------------------------------------------------

    %% 1. First derivative
    d_ecg = diff(ecg_data);
    d_ecg = [d_ecg(1); d_ecg];
    
    %% 2. Define searching window
    % x_m (maximum slope point) window (0-200 ms after T_peak)
    win_xm_start = T_peak;
    win_xm_end = min(T_peak + round(0.2*fs), length(ecg_data));
    win_xm = win_xm_start:win_xm_end;
    
    % x_r (isoelectric reference) window (200-400 ms after T_peak)
    win_xr_start = T_peak + round(0.2*fs);
    win_xr_end = min(T_peak + round(0.4*fs), length(ecg_data));
    win_xr = win_xr_start:win_xr_end;
    
    %% 3. Find x_m
    % depends on the polarity
    if ecg_data(T_peak) > 0 % positive wave
        [~, idx_min] = min(d_ecg(win_xm));
        xm = win_xm(idx_min);
        ym = ecg_data(xm);
    else                    % negative wave
        [~, idx_min] = max(d_ecg(win_xm));
        xm = win_xm(idx_min);
        ym = ecg_data(xm);
    end
    
    %% 4. Find x_r
    [value_zero, idx_zero] = min(abs(d_ecg(win_xr)));

    if value_zero == 0 
        xr = round(win_xr(idx_zero));
    else 
        xr = round(mean(win_xr)); % mean value
    end
    %xr = win_xr(idx_zero);
    yr = ecg_data(xr);
    
    %% 5. Prepare candidates (xi, yi) between x_m and x_r
    candidates = (xm+1):(xr-1);
    if isempty(candidates)
        candidates = xr;
    end
    
    %% 6. Select best option xi,yi = greatest trapezoid area
    % calculate all options
    areas = zeros(length(candidates), 1);
    
    for j = 1:length(candidates)
        xi = candidates(j);
        yi = ecg_data(xi);
        areas(j) = abs(0.5 * (ym - yi) * (2*xr - xi - xm));
    end

    % find maximum value
    [max_area, max_idx] = max(areas);
    T_end = candidates(max_idx);

    %% 7. Overview of the calculations
    t_sec = (1:length(ecg_data))/fs; %time (s)

    %info por consola
    fprintf('\n=== RESULTS ===\n');
    fprintf('T peak: %.3f s (muestra %d)\n', t_sec(T_peak), T_peak);
    fprintf('x_m: %.3f s (muestra %d)\n', t_sec(xm), xm);
    fprintf('x_r: %.3f s (muestra %d)\n', t_sec(xr), xr);
    fprintf('Detected final T: %.3f s (muestra %d)\n', t_sec(T_end), T_end);
    fprintf('Greatest area: %.2f\n', max_area);
    fprintf('RT: %.3f s (%.0f ms)\n', (T_end/fs - R_peak/fs), (T_end/fs - R_peak/fs)*1000);

    %% VISUALIZATION OF THE HOLE PROCESS
    if visual == true

        figure('Position', [100, 100, 1200, 800]);
        
        % === Subplot 1: ecg with searching windows ===
        subplot(3, 2, [1, 2]);
        
        plot(t_sec, ecg_data, 'b-', 'LineWidth', 1.5); hold on;
        
        % show T_peak
        plot(t_sec(T_peak), ecg_data(T_peak), 'r^', 'MarkerSize', 10, 'LineWidth', 1);
        
        % show searching windows
        rectangle('Position', [t_sec(win_xm(1)), min(ecg_data), ...
                              t_sec(win_xm(end))-t_sec(win_xm(1)), ...
                              max(ecg_data)-min(ecg_data)], ...
                  'EdgeColor', 'y', 'LineWidth', 1, 'LineStyle', '--', ...
                  'FaceColor', [1 1 0 0.1], 'FaceAlpha', 0.3);
        rectangle('Position', [t_sec(win_xr(1)), min(ecg_data), ...
                              t_sec(win_xr(end))-t_sec(win_xr(1)), ...
                              max(ecg_data)-min(ecg_data)], ...
                  'EdgeColor', 'm', 'LineWidth', 1, 'LineStyle', '--', ...
                  'FaceColor', [1 0 1 0.1], 'FaceAlpha', 0.3);
        
        % show x_m y x_r
        plot(t_sec(xm), ym, 'ko', 'MarkerSize', 8, 'LineWidth', 1);
        plot(t_sec(xr), yr, 'ks', 'MarkerSize', 8, 'LineWidth', 1);
        
        xlabel('Time (s)');
        ylabel('Amplitude (u)');
        title('ECG signal and searching windows')
        legend('ECG data', 'T peak', 'x_m', 'x_r', ...
               'Location', 'best');
        grid on;
        
        % === Subplot 2: derivative and x_m ===
        subplot(3, 2, 3);
        plot(t_sec, d_ecg, 'r-', 'LineWidth', 1.5); hold on;
        plot(t_sec(xm), d_ecg(xm), 'ko', 'MarkerSize', 8, 'LineWidth', 1);
        
        % reference line y=0
        yline(0, 'k--', 'LineWidth', 1);
        
        xlabel('Time (s)');
        ylabel('Derivative (u/s)');
        title('ECG signal derivative');
        legend('dECG/dt', 'Point x_m', 'Location', 'best');
        grid on;
        xlim([t_sec(xm)-0.1, t_sec(xr)+0.1]);
        
        % === Subplot 3: position vs trapezoid area ===
        subplot(3, 2, 4);
        
        plot(candidates/fs, areas, 'bo-', 'LineWidth', 1, 'MarkerSize', 2);
        hold on;
        
        % show max value
        plot(T_end/fs, max_area, 'r*', 'MarkerSize', 15, 'LineWidth', 1.5);
        
        xlabel('Candidates position (s)');
        ylabel('Area'); %ylim([0, max_area+1])
        title(sprintf('Trapezoid area (Max: %.2f en %.3f s)', max_area, T_end/fs));
        legend('areas', 'Final T point (max area)', 'Location', 'best');
        grid on;
        
        % === Subplot 4: show biggest trapezoid ===
        subplot(3, 2, [5, 6]);

        % show data between T_peak and x_r
        segment_idx = T_peak:(xr+0.2);
        plot(t_sec(segment_idx), ecg_data(segment_idx), 'b-', 'LineWidth', 2); hold on;
        
        yline(0, 'k--', 'LineWidth', 1); %y = 0

        % mark the 4 edges and the isoelectric reference point
        plot(t_sec(xm), ym, 'ro', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor', 'r');
        plot(t_sec(xr), ym, 'go', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor', 'g');
        plot(t_sec(xr), yr, 'yo', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor', 'y');
        plot(t_sec(xr), ecg_data(T_end), 'mo', 'MarkerSize', 10, 'LineWidth', 1, 'MarkerFaceColor', 'm');
        plot(t_sec(T_end), ecg_data(T_end), 'k*', 'MarkerSize', 12, 'LineWidth', 1, 'MarkerFaceColor', 'k');
        
        % legends
        text(t_sec(xm), ym, '(x_m, y_m)', ... %x_m
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
        text(t_sec(xr), ym, '(x_r, y_m)', ... 
             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10);
        text(t_sec(xr), yr, '(x_r, y_r)', ... %x_r
             'VerticalAlignment', 'top', 'FontSize', 10);
        text(t_sec(xr), ecg_data(T_end), '(x_r, y_i)', ...
             'VerticalAlignment', 'top', 'FontSize', 10);
        text(t_sec(T_end), ecg_data(T_end), sprintf('T_{end} = (x_i, y_i)\n[%.3f s]', t_sec(T_end)), ...
             'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'k');
        
        % lines between edges
        line([t_sec(xm), t_sec(xm)], [ym, yr], 'Color','black','LineStyle','--')
        line([t_sec(T_end), t_sec(T_end)], [ecg_data(T_end), yr], 'Color','black','LineStyle','--')

        line([t_sec(xr), t_sec(xr)], [yr, ecg_data(T_end)], 'Color','black','LineStyle','--')

        % draw trapezoid
        vertices_x = [t_sec(xm), t_sec(xr), t_sec(xr), t_sec(T_end)];
        vertices_y = [ym, ym, ecg_data(T_end), ecg_data(T_end)];
        
        patch(vertices_x, vertices_y, [0 1 1], 'FaceAlpha', 0.1 , 'EdgeColor', 'cyan', 'LineWidth', 1);

        % calculate and show the area
        area_value = 0.5 * (ym - ecg_data(T_end)) * (2*xr - T_end - xm);
        text(mean([t_sec(xm), t_sec(xr)]), mean([ym, yr]), ...
             sprintf('Area = %.2f', area_value), ...
             'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
        
        xlabel('Time (s)'); xlim([t_sec(T_peak)-0.15 t_sec(xr)+0.1])
        ylabel('Amplitude (u)'); %ylim([0-ecg_filtered(T_peak)*0.5 ecg_filtered(T_peak)+ecg_filtered(T_peak)*0.2])
        title(sprintf('Best option - T final point %.3f s', t_sec(T_end)));
        legend('T segment', 'x_m (max slope)', 'x_r (isoelectric reference)', 'T_{end} found', 'Area', ...
               'Location', 'west');
        grid on;
    end
    
    
end