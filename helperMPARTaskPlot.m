function scnplt = helperMPARTaskPlot(type,varargin)
% This function is only in support of MPARSearchTrackExample. It may be
% removed in a future release.


co = get(groot,'defaultAxesColorOrder');

switch type
    case 'initialize'
        % plot search grid
        figure(1);
        scanangles = varargin{1};
        azscanangles = varargin{2};
        maxrng = varargin{3};
        beamw = varargin{4};
        tgtpos = varargin{5};
        x1 = varargin{6};
        y1 = varargin{7};
        z1 = varargin{8};
        subplot(211);
        ploth.bgrid = scatter(scanangles(1,:),scanangles(2,:),20,'MarkerEdgeColor',co(1,:));
        xlabel('Azimuth Angles (deg)');
        ylabel('Elevation Angles (deg)');
        title('Search Beam Grid');
        axis equal;
        axis off;
        hold on;
        ploth.cbmarker = plot(scanangles(1,1),scanangles(2,2),'o','MarkerSize',5,...
            'MarkerEdgeColor',co(2,:),'MarkerFaceColor',co(2,:));
        xlabel('Azimuth (deg)');
        ylabel('Elevation (deg)');
        legend('Beam grid','Beam','Location','NorthWest','FontSize',6);

        % plot search area
        subplot(212);
        angmin = azscanangles(1)-beamw/2;
        angmax = azscanangles(end)+beamw/2;
        sareax = [0 maxrng*cosd([angmin azscanangles angmax]) 0];
        sareay = [0 maxrng*sind([angmin azscanangles angmax]) 0];
        ploth.bgp = patch(sareax,sareay,co(1,:),'FaceAlpha',1);
        title('Radar Azimuth Coverage');
        axis square;
        axis off;
        hold on;
        
        text(0.1*maxrng*cosd(angmin-90),0.1*maxrng*sind(angmin-90),'0','Rotation',angmin);
        text(0.8*maxrng*cosd(angmin)+0.1*maxrng*cosd(angmin-90),...
            0.8*maxrng*sind(angmin)+0.1*maxrng*sind(angmin-90),num2str(maxrng/1000),'Rotation',angmin);
        text(0.2*maxrng*cosd(angmin)+0.1*maxrng*cosd(angmin-90),...
            0.2*maxrng*sind(angmin)+0.1*maxrng*sind(angmin-90),'Range (km)','Rotation',angmin);
        
        text(1.01*maxrng*cosd(azscanangles(1)),1.01*maxrng*sind(azscanangles(1)),...
            num2str(azscanangles(1)),'Rotation',azscanangles(1));
        text(1.01*maxrng*cosd(azscanangles(end)),1.01*maxrng*sind(azscanangles(end)),...
            num2str(azscanangles(end)),'Rotation',azscanangles(end));
        

        currentaz = scanangles(1,1);
        sbeamx = [0 maxrng*cosd(currentaz+beamw/2*[-1 0 1]) 0];
        sbeamy = [0 maxrng*sind(currentaz+beamw/2*[-1 0 1]) 0];
        ploth.currentp = patch(sbeamx,sbeamy,co(2,:),'FaceAlpha',1);

        tgtmarker(1) = plot(tgtpos(1,1),tgtpos(2,1),'k^','MarkerSize',5,'MarkerFaceColor','k');
        %tgtmarker(2) = plot(tgtpos(1,2),tgtpos(2,2),'k^','MarkerSize',5,'MarkerFaceColor','k');
        tgttrack(1) = animatedline('Color','k');
        %tgttrack(2) = animatedline('Color','k');
        ploth.tgtmarker = tgtmarker;
        ploth.tgttrack = tgttrack;
        
        legend([ploth.currentp,ploth.tgtmarker(1)],'Beam','Target','Location','None',...
            'FontSize',6);
        ploth.azcovax = gca;
        ploth.azcovax.Position(4) = ploth.azcovax.Position(4)*1.1;
        hold off;
        figure(2);
        % plot target 1
        zmarkersz = 5;
        %zpossz = 200;
        %subplot(313);
        hold on;
        view(3);
        axis([35000 50000 1000 10000 -50 2000]);
        
        %ploth.tgtmarkerz(1) = plot3([tgtpos(1,1),tgtpos(1,2)],[tgtpos(2,1),tgtpos(2,2)],[tgtpos(3,1),tgtpos(3,2)],'k^','MarkerSize',zmarkersz,...
          %  'MarkerEdgeColor','k','MarkerFaceColor','k');
        
        tgtmarkerz(1) = plot3(tgtpos(1,1),tgtpos(2,1),tgtpos(3,1),'k^','MarkerSize',zmarkersz,...
            'MarkerEdgeColor','k','MarkerFaceColor','k');
        plot3(x1,y1,z1,'k.');
        %tgtmarkerz(2) = plot3(tgtpos(1,2),tgtpos(2,2),tgtpos(3,2),'k^','MarkerSize',zmarkersz,...
        %    'MarkerEdgeColor','k','MarkerFaceColor','k');
        %{
        tgtmarkerz(3) = plot3(tgtpos(1,3),tgtpos(2,3),tgtpos(3,3),'k^','MarkerSize',zmarkersz,...
            'MarkerEdgeColor','k','MarkerFaceColor','k');
        %}
         plt.tgtmarkerz = tgtmarkerz;
         
        hold off;
        title('Target 1');
        %ploth.tgtmarkerzax(1) = get(ploth.tgtmarkerz,'Parent');
        %set(ploth.tgtmarkerzax(1),'XTickLabel',get(ploth.tgtmarkerzax(1),'XTick')/1000);
        %set(ploth.tgtmarkerzax(1),'YTickLabel',get(ploth.tgtmarkerzax(1),'YTick')/1000);
        %set(ploth.tgtmarkerzax(1),'Color',co(1,:));
        xlabel('X (m)');
        ylabel('Y (m)');
        %axis([tgtpos(1,1)+zpossz*[-3 1],tgtpos(2,1)+zpossz*[-1 1]]);
        
        grid on;
        hold on;
        tgtdmarkerz(1) = plot3(nan,nan,nan,'ko','MarkerSize',zmarkersz,...
            'MarkerEdgeColor',co(2,:),'MarkerFaceColor',co(2,:));
        %tgtdmarkerz(2) = plot3(nan,nan,nan,'k^','MarkerSize',zmarkersz,...
         %   'MarkerEdgeColor',co(2,:),'MarkerFaceColor',co(2,:));
        %{
        tgtdmarkerz(3) = plot3(nan,nan,nan,'k+','MarkerSize',zmarkersz,...
          'MarkerEdgeColor',co(2,:),'MarkerFaceColor',co(2,:));
            %}
        plt.tgtdmarkerz = tgtdmarkerz;
        tgttrackz(1) = animatedline('Color',co(3,:),'MarkerSize',zmarkersz,'Marker','+');
        %tgttrackz(2) = animatedline('Color',co(3,:),'MarkerSize',zmarkersz,'Marker','s');
        %tgttrackz(3) = animatedline('Color',co(3,:),'MarkerSize',zmarkersz,'Marker','*');
        plt.tgttrackz = tgttrackz;
        %axis([tgtpos(1,1)+zpossz*[-3 1],tgtpos(2,1)+zpossz*[-1 1]]);
        legend('Target1','Trajectory1','Detection1','Track1','Location','NorthWest','FontSize',6);
        hold off;
        scnplt.ploth = ploth;
        scnplt.plt = plt;
    case 'update'
        scnplt = varargin{1};
        current_job = varargin{2};
        maxrng = varargin{3};
        beamw = varargin{4};
        tgtpos = varargin{5};
        
        scnplt.plt.tgtmarkerz(1).XData = tgtpos(1,1);
        scnplt.plt.tgtmarkerz(1).YData = tgtpos(2,1);
        scnplt.plt.tgtmarkerz(1).ZData = tgtpos(3,1);
        %ploth.tgtmarkerz(2).XData = tgtpos(1,2);
        %ploth.tgtmarkerz(2).YData = tgtpos(2,2);
        %ploth.tgtmarkerz(2).ZData = tgtpos(3,2);
        %{
        ploth.tgtmarkerz(3).XData = tgtpos(1,3);
        ploth.tgtmarkerz(3).YData = tgtpos(2,3);
        ploth.tgtmarkerz(3).ZData = tgtpos(3,3);
        %}
        scnplt.ploth.tgtmarker(1).XData = tgtpos(1,1);
        scnplt.ploth.tgtmarker(1).YData = tgtpos(2,1);
        %ploth.tgtmarker(2).XData = tgtpos(1,2);
        %ploth.tgtmarker(2).YData = tgtpos(2,2);
        currentaz = current_job.BeamDirection(1);
        currentel = current_job.BeamDirection(2);
        scnplt.ploth.cbmarker.XData = currentaz;
        scnplt.ploth.cbmarker.YData = currentel;
        sbeamx = [0 maxrng*cosd(currentaz+beamw/2*[-1 0 1]) 0];
        sbeamy = [0 maxrng*sind(currentaz+beamw/2*[-1 0 1]) 0];
        scnplt.ploth.currentp.XData = sbeamx;
        scnplt.ploth.currentp.YData = sbeamy;
        
        switch current_job.JobType
            case 'Search'
                scnplt.ploth.cbmarker.MarkerEdgeColor = co(2,:);
                scnplt.ploth.cbmarker.MarkerFaceColor = co(2,:);
                scnplt.ploth.currentp.FaceColor = co(2,:);
            case 'Confirm'
                scnplt.ploth.cbmarker.MarkerEdgeColor = co(3,:);
                scnplt.ploth.cbmarker.MarkerFaceColor = co(3,:);
                scnplt.ploth.currentp.FaceColor = co(3,:);
                trackid = current_job.TrackID;
                allTracks = varargin{6};
                dm = varargin{7};
                xt = allTracks(trackid).State(1);
                yt = allTracks(trackid).State(3);
                zt = allTracks(trackid).State(5);
                addpoints(scnplt.plt.tgttrackz(trackid),xt,yt,zt);
                [xd,yd,zd] = sph2cart(deg2rad(dm(1)),deg2rad(dm(2)),dm(3));
                scnplt.plt.tgtdmarkerz(trackid).XData = xd;
                scnplt.plt.tgtdmarkerz(trackid).YData = yd;
                scnplt.plt.tgtdmarkerz(trackid).ZData = zd;
                %addpoints(ploth.tgttrackz(trackid),xt,yt);

            case 'Track'
                scnplt.ploth.cbmarker.MarkerEdgeColor = co(4,:);
                scnplt.ploth.cbmarker.MarkerFaceColor = co(4,:);
                scnplt.ploth.currentp.FaceColor = co(4,:);
                trackid = current_job.TrackID;
                allTracks = varargin{6};
                dm = varargin{7};
                xt = allTracks(trackid).State(1);
                yt = allTracks(trackid).State(3);
                zt = allTracks(trackid).State(5);
                addpoints(scnplt.plt.tgttrackz(trackid),xt,yt,zt);
                [xd,yd,zd] = sph2cart(deg2rad(dm(1)),deg2rad(dm(2)),dm(3));
                %ploth.tgtmarkerz(trackid).XData = tgtpos(1,trackid);
                %ploth.tgtmarkerz(trackid).YData = tgtpos(2,trackid);
                scnplt.plt.tgtdmarkerz(trackid).XData = xd;
                scnplt.plt.tgtdmarkerz(trackid).YData = yd;
                scnplt.plt.tgtdmarkerz(trackid).ZData = zd;
                %addpoints(ploth.tgttrackz(trackid),xt,yt);
        end
end

drawnow;




















































%Pranava K Bhat 1BY17TE026
