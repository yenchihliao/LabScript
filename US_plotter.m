% Script for producing plots based on DeepEdge output

clear;

%%%%%%%%%%%%% Specify settings  here %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
config = struct();

% Number of frames for animated plots

config.frames = 30;

% Millimeters per pixel: use the "Measure 1 cm" command under the
% "commands" dropdown menu in DeepEdge

config.mpp = 0.166; %Tml
%config.mpp = 0.173; %Nnw
%config.mpp = 0.151; %Ksv

% xlims and ylims

config.xlims = [60 160];
config.ylims = [-100 0];

config.ema_xlims = [-55 5];
config.ema_ylims = [-30 30];

%%%%%%%%%%%%%%%%%%%%%%%End of user settings%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select working directory
parent_directory = uigetdir('',"Select working directory");
config.parent_directory = parent_directory;

cd(parent_directory);

if ~isfolder('fig')
    mkdir('fig');
end
figure_directory = fullfile(parent_directory,'fig');

% Isolate contour files
contour_files = dir('*.mat');
contour_file_names = {contour_files.name};


% Choose what to plot
%{ WIP
    '2D plot n points (single trial)', ...
    '2D plot n points (averaged)', ...
    '2D plot Animated (averaged)',
%}
all_plot_types = {
    '2D plot (single trial)', ...
    '2D plot (all trials)', ...
    '2D plot (averaged)', ...
    '2D plot Animated (single trial)', ...
    'Compare Averages', ...
    'Align EMA', ...
    'Diff'
    };
[plot_index,plot_tf] = listdlg('PromptString','Select plot types.','ListString',all_plot_types);

plot_types = all_plot_types(plot_index);

if any(ismember({'2D plot (single trial)' '2D plot (all trials)' '2D plot (averaged)' 'Compare Averages' 'Align EMA'},plot_types))

    t = str2double(inputdlg('Enter the desired time for the plots (0 ~ 1):'));
    while isnan(t) || t < 0 || t > 1
        t = str2double(inputdlg('Error - you must enter a number between 0.0 and 1.0. Try again:'));
    end
end

for m = 1:length(plot_types)
    switch plot_types{m}
        % Select multiple distinct file(s) at the same time, plotted separatly
        case '2D plot (single trial)'
            [contour_indexes,~] = listdlg(...
                                    'PromptString','Select which files for 2D plot (single trial).',...
                                    'ListString',contour_file_names);
            target_contours = contour_file_names(contour_indexes);
            for c_idx = 1:length(target_contours)
                c_file = target_contours{c_idx};
                save_file_name = c_file(1:end-4); % default save file name is just the same as the TSV file name
                c_data = importdata(c_file);

                cd(figure_directory);
                Plot2DToFile(c_data,save_file_name,t,config);

                cd(parent_directory);
            end
        % Selcet multiple distinct file at the same time, plotted together
        case '2D plot (all trials)'
            [contour_indexes,~] = listdlg('PromptString','Select which files for 2D plot together.','ListString',contour_file_names);
            target_contours = contour_file_names(contour_indexes);
            PlotAll2File(target_contours, t, config)

        case '2D plot (averaged)'
            [contour_indexes,~] = listdlg('PromptString','Select which files for 2D plot (averaged).','ListString',contour_file_names);
            target_contours = contour_file_names(contour_indexes);
            all_data = {};
            for c_idx = 1:length(target_contours)
                c_file = target_contours{c_idx};
                save_file_name = c_file(1:end-4); % default save file name is just the same as the TSV file name
                c_data = importdata(c_file);
                all_data{c_idx} = c_data;
            end
            cd(figure_directory);
            PlotAverageContour(all_data,t,config);
            cd(parent_directory);
        % case '2D plot n points (single trial)'
        % case '2D plot n points (averaged)'
        case '2D plot Animated (single trial)'
            [contour_indexes,~] = listdlg('PromptString','Select which files for 2D Animated (single trial).','ListString',contour_file_names);
            target_contours = contour_file_names(contour_indexes);
            for c_idx = 1:length(target_contours)
                c_file = target_contours{c_idx};
                save_file_name = c_file(1:end-4); % default save file name is just the same as the TSV file name
                c_data = importdata(c_file);

                cd(figure_directory);
                Animate2DIndividual(c_data,save_file_name, config)
                cd(parent_directory);
            end
        % case '2D plot Animated (averaged)'
        case 'Compare Averages'
            [contour_indexes,~] = listdlg('PromptString','Select first set to average.','ListString',contour_file_names);
            target_contoursA = contour_file_names(contour_indexes);
            [contour_indexes,~] = listdlg('PromptString','Select second set to average.','ListString',contour_file_names);
            target_contoursB = contour_file_names(contour_indexes);
            all_dataA = {}; all_data_B = {};
            for c_idx = 1:length(target_contoursA)
                c_file = target_contoursA{c_idx};
                save_file_name = c_file(1:end-4); % default save file name is just the same as the TSV file name
                c_data = importdata(c_file);
                all_dataA{c_idx} = c_data;
            end

            for c_idx = 1:length(target_contoursB)
                c_file = target_contoursB{c_idx};
                save_file_name = c_file(1:end-4); % default save file name is just the same as the TSV file name
                c_data = importdata(c_file);
                all_dataB{c_idx} = c_data;
            end
            cd(figure_directory);
            PlotAverageContoursComparison(all_dataA,all_dataB,t,config);
            cd(parent_directory);
        case 'Align EMA'
            [contour_indexes,~] = listdlg('PromptString','Select which files for EMA comparison (averaged).','ListString',contour_file_names);
            target_contours = contour_file_names(contour_indexes);
            all_data = {};

            [fName,path] = uigetfile('*.mat',"Select EMA average file");

            EMA_ave = importdata(fullfile(path,fName));

            for c_idx = 1:length(target_contours)
                c_file = target_contours{c_idx};
                save_file_name = c_file(1:end-4); % default save file name is just the same as the TSV file name
                c_data = importdata(c_file);
                all_data{c_idx} = c_data;
            end
            cd(figure_directory);
            AlignedEMAPlot(all_data, EMA_ave, t, config);
            cd(parent_directory);
        % Select multiple file from any file at any time, plotted together
        case 'Diff'
            % To store user inputs
            contours = {};
            times = [];

            selected = '';
            t = '0.5';
            while true
                name = sprintf('Selected %d files', length(contours));
                [selected,~] = listdlg(...
                                        'Name', name,...
                                        'PromptString','Select a file or cancel',...
                                        'SelectionMode','single',...
                                        'ListString',[contour_file_names]);
                % If user cancels, break the loop
                if isempty(selected)
                    break
                end
                contour = contour_file_names(selected);
                contours = [contours, contour];
                t = str2double(...
                    inputdlg(...
                        'Enter the desired time for the plots (0 ~ 1):'),...
                        'Select time', 1,...
                        {num2str(t)}... % Use the previous time as default
                    );
                while isnan(t) || t < 0 || t > 1
                    t = str2double(inputdlg('Error - you must enter a number between 0.0 and 1.0. Try again:'));
                end
                times = [times, t];
            end
            PlotAll2File(contours, times, config)
    end
    cd(parent_directory)
end




%%%%%
function [contourX, contourY] = GetContour(data,time,config)
    mpp = config.mpp;

    n_frames = length(data);
    t_frame = round(time*n_frames);
    if t_frame == 0
        t_frame = 1;
    end
    contourX = data(t_frame).XY(:,1)*mpp; % in mm
    contourY = data(t_frame).XY(:,2)*-1*mpp; % in mm (inverted because y scale is inverted)

end

function [ave_contour,all_contours] = GetAverageContour(all_data,time,config)
    mpp = config.mpp;

    n_contours = length(all_data);

    all_contours = [];
    for c_idx = 1:n_contours
        c_data = all_data{c_idx};
        n_frames = length(c_data);
        t_frame = round(time*n_frames);
        if t_frame == 0
            t_frame = 1;
        end
        contourX = c_data(t_frame).XY(:,1)*mpp; % in mm
        contourY = c_data(t_frame).XY(:,2)*-1*mpp; % in mm (inverted because y scale is inverted)
        all_contours = [all_contours contourX contourY];
    end

    ave_contour = [];
    for r_idx = 1:length(all_contours)
        row = all_contours(r_idx,:);
        xave = mean(row(1:2:length(row)));
        yave = mean(row(2:2:length(row)));
        ave_contour = [ave_contour; xave yave];
    end
end

%{
INPUT:
    contours: Cell array of contour names
    times: Time for all the contours. Or a List of times over each contour
    config
%}
function PlotAll2File(contours, times, config)
    fig2d = figure('Visible','off','position',[0, 0, 500,500]);
    len = length(contours);
    if length(times) ~= len
        if length(times) ~= 1
            error("PlotAll2File: Contours and times length mismatch")
        else
            % When times is a single number
            times = times * ones(1, len);
        end
    end

    for i = 1:len
        [cX, cY] = GetContour(importdata(contours{i}),times(i),config);
        display = split(contours{i}, '_');
        % e.g. ZAX Halrup 0.3
        display = strjoin([display{1}, ' ', display{2} ,' ', string(times(i))])
        plot(cX,cY, 'DisplayName', display);
        hold on
    end

    xlim(config.xlims);
    ylim(config.ylims);
    legend('Interpreter', 'none')

    cd("./fig");
    output_file_name = inputdlg('Save file as... (no extension)');
    saveas(fig2d,append(output_file_name,".png"));
    % cla;
    close(fig2d);
end

function PlotAverageContour(all_data,time,config)
    xlims = config.xlims;
    ylims = config.ylims;

    [ave_contour, all_contours] = GetAverageContour(all_data,time,config);
    fig2d = figure('Visible','off','position',[0, 0, 500,500]);

    % plot average
    avex = ave_contour(:,1);
    avey = ave_contour(:,2);

    plot(avex,avey,"b");
    hold on
    xlim(xlims);
    ylim(ylims);
    % xlabel('x-position (mm)');
    % ylabel('z-position (mm)');
    for c_idx = 1:length(all_contours(1,:))/2
        c_x = all_contours(:,1+(c_idx-1)*2);
        c_y = all_contours(:,2+(c_idx-1)*2);
        plot(c_x,c_y,'Color',[0,0,1,0.3]);
    end

    output_file_name = inputdlg('Save file as... (no extension)');
    saveas(fig2d,append(output_file_name,".png"));
    cla;
    close(fig2d);
end

function PlotAverageContoursComparison(all_dataA,all_dataB,time,config)
    xlims = config.xlims;
    ylims = config.ylims;

    [ave_contourA, all_contoursA] = GetAverageContour(all_dataA,time,config);
    [ave_contourB, all_contoursB] = GetAverageContour(all_dataB,time,config);
    fig2d = figure('Visible','off','position',[0, 0, 500,500]);

    % set 1
    % plot average
    avex = ave_contourA(:,1);
    avey = ave_contourA(:,2);

    plot(avex,avey,"b","LineWidth",1.5);
    hold on
    xlim(xlims);
    ylim(ylims);
    % xlabel('x-position (mm)');
    % ylabel('z-position (mm)');
    for c_idx = 1:length(all_contoursA(1,:))/2
        c_x = all_contoursA(:,1+(c_idx-1)*2);
        c_y = all_contoursA(:,2+(c_idx-1)*2);
        plot(c_x,c_y,'Color',[0,0,1,0.3]);
    end

    % set 2
    % plot average
    avex = ave_contourB(:,1);
    avey = ave_contourB(:,2);

    plot(avex,avey,"r","LineWidth",1.5);

    for c_idx = 1:length(all_contoursB(1,:))/2
        c_x = all_contoursB(:,1+(c_idx-1)*2);
        c_y = all_contoursB(:,2+(c_idx-1)*2);
        plot(c_x,c_y,'Color',[1,0,0,0.3]);
    end

    output_file_name = inputdlg('Save file as... (no extension)');
    saveas(fig2d,append(output_file_name,".png"));
    cla;
    close(fig2d);
end

function Plot2DToFile(data,file_name,time,config)
    fprintf("Plotting 2D individual contour for %s... \n",file_name);
    fig = Plot2DIndividual(data,time,config);
    output_file_name = sprintf("%s_%d.png",file_name,time*100);
    saveas(fig,output_file_name);
    cla;
    close(fig);

end

function fig2d = Plot2DIndividual(data,time,config)
    xlims = config.xlims;
    ylims = config.ylims;

    [cX, cY] = GetContour(data,time,config);

    fig2d = figure('Visible','off','position',[0, 0, 500,500]);
    plot(cX,cY);
    xlim(xlims);
    ylim(ylims);
    % xlabel('x-position (mm)');
    % ylabel('z-position (mm)');

    % Remove coordinates
    ax = gca;
    ax.XTick = []; % Remove numerical x-ticks
    ax.YTick = []; % Remove numerical y-ticks

    % Set axis with tongue indicator
    text(xlims(2), ylims(1), 'TT', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(xlims(1), ylims(1), 'TD', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    text(xlims(1), ylims(2), 'Upper', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(xlims(1), ylims(1), 'Lower', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

function fig2d = Plot2DAll(all_data,time,config)
    xlims = config.xlims;
    ylims = config.ylims;
    fig2d = figure('Visible','off','position',[0, 0, 500,500]);

    for i = 1:length(all_data)
        data = all_data{i};
        [cX, cY] = GetContour(data,time,config);

        plot(cX,cY,'Color','b');
        hold on
    end
    xlim(xlims);
    ylim(ylims);
    % xlabel('x-position (mm)');
    % ylabel('z-position (mm)');
end

function Animate2DIndividual(data,file_name, config)
    frames = config.frames;

    gifFile = sprintf('%s_2D_animated.gif',file_name);

    times = linspace(1/frames,1,frames);

    fprintf("Plotting %s...\n",gifFile);

    for ct_idx = 1:frames
        if ct_idx == 1
            fprintf("Animating... Frame %d of %d",ct_idx,frames)
        else
            if ct_idx > 10 && frames >= 10
                fprintf("\b\b\b\b\b\b\b\b%d of %d",ct_idx,frames)
            elseif ct_idx > 10 && frames < 10
                fprintf("\b\b\b\b\b\b\b%d of %d",ct_idx,frames)
            elseif ct_idx <= 10 && frames >= 10
                fprintf("\b\b\b\b\b\b\b%d of %d",ct_idx,frames)
            else
                fprintf("\b\b\b\b\b\b%d of %d",ct_idx,frames)
            end
        end
        c_t = times(ct_idx);
        fig = Plot2DIndividual(data,c_t,config);
        if ct_idx == 1
            exportgraphics(fig,gifFile,Append=false);
        else
            exportgraphics(fig,gifFile,Append=true);
        end

    end
end

function AlignedEMAPlot(all_data, EMA_ave, time, config)
    xlims = config.ema_xlims;
    ylims = config.ema_ylims;

    [ave_contour,~] = GetAverageContour(all_data,time,config);

    TTx = EMA_ave.TT(1); TTz = EMA_ave.TT(3);
    TBx = EMA_ave.TB(1); TBz = EMA_ave.TB(3);
    TDx = EMA_ave.TD(1); TDz = EMA_ave.TD(3);

    palLine = EMA_ave.PAL;

    palX = palLine.SIGNAL(:,1);
    palZ = palLine.SIGNAL(:,3);

    TT_rel = [TTx-TBx TTz-TBz];
    TT_dist = sqrt(TT_rel(1)^2+TT_rel(2)^2);
    TT_angle = atan2(TT_rel(2),TT_rel(1));
    TD_rel = [TDx-TBx TDz-TBz];
    TD_dist = sqrt(TD_rel(1)^2+TD_rel(2)^2);
    TD_angle = atan2(TD_rel(2),TD_rel(1));

    TT_TD_angle_diff = TT_angle - TD_angle;

    error = 20;
    for i = 2:length(ave_contour)
        temp_TB = ave_contour(i,:);
        offset_x = TBx-temp_TB(1); offset_z = TBz-temp_TB(2);
        for j = i:length(ave_contour)
            temp_TT = ave_contour(j,:);
            temp_TT_rel = [temp_TT(1)-temp_TB(1) temp_TT(2)-temp_TB(2)];
            temp_TT_dist = sqrt(temp_TT_rel(1)^2+temp_TT_rel(2)^2)-TT_dist;
            if temp_TT_dist > 5
                continue
            end

            temp_TT_angle = atan2(temp_TT_rel(2),temp_TT_rel(1));
            exp_TD_angle = temp_TT_angle + TT_TD_angle_diff;

            exp_TD = [temp_TB(1)+cos(exp_TD_angle)*TD_dist temp_TB(2)+sin(exp_TD_angle)*TD_dist];

            for k = 1:i
                temp_TD = ave_contour(k,:);

                temp_TD_dist = sqrt((temp_TD(1)-exp_TD(1))^2 + (temp_TD(2)-exp_TD(2))^2);
                if temp_TD_dist > 10
                    continue
                end
                temp_error = temp_TT_dist + temp_TD_dist;
                if temp_error <= error
                    error = temp_error;
                    offset = [temp_TB(1)-TBx temp_TB(2)-TBz];
                    rotation = temp_TT_angle - TT_angle;
                    c_idx = i;
                end
            end
        end
    end
    if error == 20
        fprintf("Can't find a good fit. \n");
        return
    end

    fig2d = figure('Visible','off','position',[0, 0, 500,500]);
    plot(palX,palZ,'k','LineWidth',1.5);
    hold on
    xlim(xlims);
    ylim(ylims);
    % xlabel('x-position (mm)');
    % ylabel('z-position (mm)');

    midsagittals_x = [TDx, TBx, TTx];
    midsagittals_z = [TDz, TBz, TTz];
    plot(midsagittals_x,midsagittals_z,'.','Color','b','MarkerSize',15);

    adjusted_curve = AdjustCurve(ave_contour,c_idx,offset,rotation);
    plot(adjusted_curve(:,1),adjusted_curve(:,2),'Color',[0,1,0,.85]);
    hold off

    output_file_name = inputdlg('Save file as... (no extension)');
    saveas(fig2d,append(output_file_name,".png"));
    cla;
    close(fig2d);
end

function adjusted_curve = AdjustCurve(curve,c_idx,offset,rotation)
    center = curve(c_idx,:);
    center_x = center(1);
    center_z = center(2);
    adjusted_curve = [];

    offset_x = offset(1);
    offset_z = offset(2);

    for i = 1:length(curve)
        if i == c_idx
            adjusted_curve = [adjusted_curve; [center_x-offset_x center_z-offset_z]];
        end
        c_x = curve(i,1); c_z = curve(i,2);
        c_dist = sqrt((c_x-center_x)^2+(c_z-center_z)^2);
        c_angle = atan2(c_z-center_z,c_x-center_x);
        adj_angle = c_angle + rotation;
        adj_pos = [center_x + cos(adj_angle)*c_dist - offset_x center_z + sin(adj_angle)*c_dist - offset_z];
        adjusted_curve= [adjusted_curve; adj_pos];
    end
end