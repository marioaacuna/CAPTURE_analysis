% this is a test script to plot and make videos of clusters from the 
% unsupervised clustering,  and body poses
% taken from tsne 
% this script needs to be run after the RUN_CAPTURE_mario.m script
%% Plot Clusters
figure(1)
hold on

events = [T.frameStart, T.frameStop];
for  iv = 1:n_events
    these_events = events(iv,1):events(iv,2);
    frs = round(these_events/10);
    plot(zvals(frs,1),zvals(frs,2), 'or', 'MarkerFaceColor', 'r')

end
hold off
%% Plot clusters

h1=figure( 'Position',[10 10 2000 1200]);
CL = cell(length(seq_c_idx),1);

ncols = 1;
nrows = 1;

for seq_ic = 1:numel(seq_cls)

    this_cls = seq_cls(seq_ic);    fprintf('ic = %i - ', this_cls)
    if make_video
        writerObj = VideoWriter(fullfile('D:/_test_label3D/videos/clusters',['cluster_', num2str(this_cls),'.avi'])); %#ok<TNMLP>
        writerObj.Quality = 50;
        writerObj.FrameRate = 30;
        open(writerObj);
    end
    these_frames = find(hierarchystruct.clustered_behavior{1}==this_cls);
    CL(seq_ic) =  {find(hierarchystruct.clustered_behavior{1}==this_cls)};
    ax = subplot(nrows, ncols, 1);
    if this_cls==0,  fprintf('\n'),continue, end

    %% the default case
    [maxval] = max(mocapstruct.markers_preproc.SpineM,[],1);
    [minval] =   min(mocapstruct.markers_preproc.SpineM,[],1);
    buffactor_axis = 1.3; %1.1
    buffactor_arena = 1.02;


    th = 0:pi/100:2*pi;
    xcent = (maxval(1)+minval(1))./2;
    ycent = (maxval(1)+minval(1))./2;

    frame_last = 0;

    marker_plot = ones(1,numel(mocapstruct.markernames));


    %% initialize the figure
    set(h1,'Color','k')

    xx = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(1,1));
    yy = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(1,2));
    zz = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{1})(1,3));
    handle_base = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);


    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid off;
    axis off
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);
    %ori:
    % zlim([200 350])
    % xlim([-340 340])
    % ylim([-340 340])

    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
    view([-22, 12]);




    xunit = 304*buffactor_arena * cos(th) + xcent;
    yunit = 304*buffactor_arena * sin(th) + ycent;
    zunit = zeros(1,numel(xunit));
    plot3(xunit,yunit,zunit,'w','linewidth',3);

%     zlim(zlimvals)
%     xlim(xlimvals)
%     ylim(ylimvals)

    zlim([-110 170])
    xlim([-140 140])
    ylim([-140 140])





    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])

    axis tight



    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);


    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width 1.1*ax_height];
    view([162, 7]);

    for lk = reshape(these_frames,1,[])%1:10:10000
        cla;
        %     fprintf('frame: %f \n', lk)

        plot3(xunit,yunit,zunit,'w','linewidth',3);

        ind_to_plot = lk;

        %% Plot markers that are tracked in the frame
        set(gca,'Nextplot','ReplaceChildren');
        handles_here = cell(1,numel(mocapstruct.markernames));
%         title(ax, str_title, 'Color','w')
        for jj = 1:numel(mocapstruct.markernames)
            % don't plot markers that drop out
            if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
                if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
                    xx = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                    yy = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                    zz = squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                    handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',5);



                    hold on
                    marker_plot(jj) = 1;
                else
                    marker_plot(jj) = 0;
                end

            end
        end

        %% plot the links between markers
        for mm = 1:numel(mocapstruct.links)
            if numel(mocapstruct.links{mm})
                if (ismember(mocapstruct.links{mm}(1),1:numel(mocapstruct.markernames)) && ismember(mocapstruct.links{mm}(2),1:numel(mocapstruct.markernames)))
                    if (marker_plot(mocapstruct.links{mm}(1)) == 1 && marker_plot(mocapstruct.links{mm}(2)) == 1)

                        xx = [squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                            squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                        yy = [squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                            squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                        zz = [squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                            squeeze(mocapstruct.markers_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
                        line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',1);
                    end

                end
            end
        end

        %new
        zlim([-20 80])
        xlim([-100 100])
        ylim([-120 120])



        %     zlim(zlimvals)
        %     xlim(xlimvals)
        %     ylim(ylimvals)
        set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])

        axis off

        drawnow
        hold off

        sgtitle(['Cluster : ', num2str(this_cls)], 'color','w')
        frame_last = lk;
        Fi= getframe(h1);
        writeVideo(writerObj, Fi)



        try
            text(-80.1,-50.10,num2str(frame_last), 'color', 'w')
        catch
            keyboard
        end
        %         M(find(frame_inds == lk)) =  getframe(gcf);

        %clf
    end
    %     close(writerObj);
    
    if make_video
        close(writerObj);
    end



end

%% Plot pose estimations

make_video = 0;
h2 = figure( 'Position',[10 10 2000 1200]);

if make_video
    writerObj = VideoWriter(fullfile('D:/_test_label3D/videos/poses',['Poses_', '.avi']));
    writerObj.Quality = 50;
    writerObj.FrameRate = 100;
    open(writerObj);
end

for ic = 1:numel(cls)



    this_cls = cls(ic);
    fprintf('ic = %i - ', this_cls)
    frame_inds = find(analysisstruct.annot_reordered{end}==this_cls);

    frame_last = 0;

    marker_plot = ones(1,numel(mocapstruct.markernames));


    %% initialize the figure
    set(h2,'Color','k')

    xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,1));
    yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,2));
    zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{1})(1,3));
    handle_base = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{1},'MarkerFaceColor',mocapstruct.markercolor{1},'MarkerSize',6);


    ax = gca;
    axis(ax,'manual')
    set(gca,'Color','k')
    grid on;
    set(gca,'Xcolor',[1 1 1 ]);
    set(gca,'Ycolor',[1 1 1]);
    set(gca,'Zcolor',[1 1 1]);

    zlim([-110 170])
    xlim([-140 140])
    ylim([-140 140])


    %     zlim([-210 270])
    %     xlim([-240 240])
    %     ylim([-240 240])
    %

    set(gca,'XTickLabels',[],'YTickLabels',[],'ZTickLabels',[])
    view([-22, 12]);
    fprintf('n frames %i \n',numel(frame_inds));

    for lk = reshape(frame_inds,1,[])%1:10:10000

        cla;


        %mocapstruct.links{20} = [];
        %       mocapstruct.links{22} = [];

        ind_to_plot = lk;

        %% Plot markers that are tracked in the frame
        set(gca,'Nextplot','ReplaceChildren');
        handles_here = cell(1,numel(mocapstruct.markernames));
        %         title(texthere, 'Color','w')
        for jj = 1:numel(mocapstruct.markernames)
            % don't plot markers that drop out
            if ~isnan(sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2))
                if (~sum(mocapstruct.markers_preproc.(mocapstruct.markernames{jj})(ind_to_plot,:),2) == 0)
                    xx = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,1));
                    yy = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,2));
                    zz = squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{jj})(ind_to_plot,3));
                    handles_here{jj} = line(xx,yy,zz,'Marker','o','Color',mocapstruct.markercolor{jj},'MarkerFaceColor',mocapstruct.markercolor{jj},'MarkerSize',9);



                    hold on
                    marker_plot(jj) = 1;
                else
                    marker_plot(jj) = 0;
                end

            end
        end

        %% plot the links between markers
        for mm = 1:numel(mocapstruct.links)
            if numel(mocapstruct.links{mm})
                if (ismember(mocapstruct.links{mm}(1),1:numel(mocapstruct.markernames)) && ismember(mocapstruct.links{mm}(2),1:numel(mocapstruct.markernames)))
                    if (marker_plot(mocapstruct.links{mm}(1)) == 1 && marker_plot(mocapstruct.links{mm}(2)) == 1)

                        xx = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,1)) ...
                            squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,1)) ];
                        yy = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,2)) ...
                            squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,2))];
                        zz = [squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(1)})(ind_to_plot,3)) ...
                            squeeze(mocapstruct.markers_aligned_preproc.(mocapstruct.markernames{mocapstruct.links{mm}(2)})(ind_to_plot,3))];
                        line(xx,yy,zz,'Color',mocapstruct.markercolor{mocapstruct.links{mm}(1)},'LineWidth',3);
                    end

                end
            end
        end

        drawnow

        hold off

        sgtitle(['Cluster : ', num2str(this_cls)], 'color','w')
        frame_last = lk;
        if make_video
            Fi= getframe(h2);
            writeVideo(writerObj, Fi)
        end

        %         frame_last = lk;

        %         M(find(frame_inds == lk)) =  getframe(gcf);

        %clf
    end


end

if make_video
    close(writerObj);
end
