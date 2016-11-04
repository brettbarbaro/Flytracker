function flytracker

% this function will run all video files in the current directory and
% return the results in a single csv file for import into Excel

    % MASTER VARIABLES

    version=5.1 %MAKE SURE TO UPDATE THIS EVERY TIME A NEW VERSION IS MADE!

                %5.1 - addressing bug where the program doesn't work for 4
                %vials.

                %NOTE - make sure to let some flies get all the way to the top
                %for each segment. Otherwise, the box-drawing algorithm won't
                %work properly. It's best to have a positive control vial and
                %to let each segment go to at least 4 seconds. 5 probably
                %better. Note also that after the flies get to the top, they
                %won't be easy to detect due to the background subtraction
                %method using the last frame.

    fps=30;  %this must match the frames per second of the input video

    sens = .4;  %sensitivity of thresholding for segmentation - .4 for most
                %but might want to adjust for very dim or bright videos.
                %doesn't affect thresholding for fly finding, because that is
                %done later with graythresh on background-subtracted images.
                %This was introduced when analyzing Marianne's 320fps videos,
                %which were very dark due to being high-speed.

    framerate=fps/6; %NOTE: framerate must divide evenly into fps. 
                % Taking a frame every 1/6 of a second works pretty well for flies.

    step=framerate*3;

fprintf('\nWelcome to Flytracker!\n\n');
fprintf('Flytracker will analyze all .avi and .mov files in the current folder.\n\n');
fprintf('For best results, vials should be evenly spaced\nwith healthy control flies in the outside vials.\n\n');

numvials = input('How many vials are being compared? ');

verbose = input('Do you want to see the details? (0=no, 1=yes) ');

% if size(varargin,2)==0
%     numvials=6;
%     fprintf('\ndefault number of vials=6\n\n')
% end
% 
% if size(varargin,2)==1
%     numvials=varargin{1};
%     fprintf('\nnumber of vials=%d\n\n',numvials)
% end

currentDirectory = pwd; 
[~, deepestFolder] = fileparts(currentDirectory);
filename=strcat(deepestFolder,'_OUTPUT.csv');


output=[];

avi=dir('*.AVI');
mov=dir('*.mov');
dirlist = [avi; mov];

for i=1:numel(dirlist)
    
        tic

    %BEGIN SEGMENT

    fprintf('\nreading ');
    fprintf(dirlist(i).name);
    fprintf('\n');
    a=VideoReader(dirlist(i).name);


    %first, output the first frame, which should have the labeled vials

    h=imshow(read(a,1));
    saveas(h,strcat(dirlist(i).name,'_v',num2str(version),'_FIRSTFRAME.jpg'));


    % replace bad characters in video name with acceptable chars
    % excel won't open files with / or : in them
    dirlist(i).name = strrep(dirlist(i).name,'/','_');
    dirlist(i).name = strrep(dirlist(i).name,':','_');  


    diff=zeros(fix(a.NumberOfFrames/step),2);

    %take a sample of the video, every <step> frame
    %and compare adjacent samples by subtraction
    %sum up the differences between pixels
    %this gives overview of video
    %experimental sequences start with low number runs

    y=0;

    disp('creating overview')
    fprintf('\nat frame:');
    fprintf('          \n');

    for x=1:step:a.NumberOfFrames-1-step
        if mod(x,(step*100))==1
            fprintf(1,'\b\b\b\b\b\b\b\b\b\b%10.0f',x);
        end
        y=y+1;
        diff(y,1)=x;
        c=im2bw(read(a,x));
        d=im2bw(read(a,x+step));
        e=c-d;
        f=abs(e);
        g=(sum(sum(f)));
        diff(y,2)=g;
    end

    %the median value of the differences should be approximately
    %the difference between frames of fly-counting runs, provided that the
    %fly-counting frames outnumber the in-between frames
    %this should probably be improved in later versions to extract the relevant
    %data more precisely

    %double the median number is taken as the threshhold for in-between frames,
    %and half the median is taken as the threshhold for aberrant in-between
    %frames (e.g. when there's nothing there at all)

    diffs=diff(:,2);
    index=median(diffs);
    doubleindex=2*index;
    halfindex=0;

    diffs(size(diffs,1)-1)=3*doubleindex; %creating an endpoint so edge detector below will find it

    y=1;
    k=0;
    roughseg=0;

    for x=3:size(diffs)-2
        if (k==0 && diffs(x-1)>doubleindex && diffs(x)>halfindex && diffs(x+1)>halfindex && diffs(x+2)>halfindex && diffs(x)<doubleindex && diffs(x+1)<doubleindex && diffs(x+2)<doubleindex)
           roughseg(y,1)=1+(x-1)*step;
           k=1;
        end
        if (k==1 && diffs(x+1)>doubleindex && diffs(x)>halfindex && diffs(x-1)>halfindex && diffs(x-2)>halfindex && diffs(x)<doubleindex && diffs(x-1)<doubleindex && diffs(x-2)<doubleindex)
           roughseg(y,2)=1+(x-1)*step;  %#ok<*AGROW>
           y=y+1;
           k=0;
        end
    end

    seg=roughseg;
    maxseglength=0; % maxseglength is the length of the longest segment, in 
                    % seconds, and is used to create the proper size matrices
                    % below. not sure if it's really necessary - maybe should
                    % allocate matrices of predetermined length, like 30. in
                    % its current form, this program outputs a matrix of length
                    % 30 anyway.

    for x=1:size(roughseg,1)
        y=roughseg(x,1);
        g=index;
        while g<doubleindex && y>1
            c=im2bw(read(a,y),sens);
            d=im2bw(read(a,y-1),sens);
            e=c-d;
            f=abs(e);
            g=(sum(sum(f)));
            y=y-1;
        end
        seg(x,1)=y;
        seglength=fix((seg(x,2)-y)/fps);
        if seglength>maxseglength
            maxseglength=seglength;
        end
    end

        % "seg" is a matrix with a row for each segment of the video. 
        % the first column is the start point, the second is the end.

    % USE THIS CODE IF YOU WANT TO OUTPUT A VIDEO
    % 
    % for x=1:size(seg,1)
    %     segment=strcat(video,'_v',num2str(version),'_seg',num2str(x),'.avi')
    %     aviobj=avifile(segment);
    %     for y=seg(x,1):seg(x,2)
    %         F=read(a,y);
    %         aviobj=addframe(aviobj,F);
    %     end
    %     aviobj=close(aviobj);
    % end

    %END SEGMENT


    segments=size(seg,1);
    ahpsperseg=zeros(maxseglength,0); % for AVERAGE HEIGHTS
    thpsperseg=zeros(maxseglength,0); % for TOTAL HEIGHTS

    totsecs=0;
    
    for x=1:maxseglength
        totsecs(x,1)=x;
    end



    %START ANALYZE

    fprintf('\n\nanalyzing segments\n');
    
    for segment=1:segments

        fprintf('\nsegment %d\n\n',segment)

        segbeg=seg(segment,1);
        segend=seg(segment,2);
        % define the beginning and ending frames for the segment being analyzed

        if verbose==1
            h=imshow(read(a,segbeg));
            saveas(h,strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_frame1.jpg'));
        end

        % FIND BACKGROUND
        % takes a VideoReader object, takes its inverse (so that flies are white
        % instead of black - black has a numerical value of zero, which makes
        % it harder to work with), sums up every <framerate>th frame, and then
        % calculates their average. This should yield the elements of the
        % segment that are static, so that they can be subtracted out. 
        % NOTE: IF A FLY IS
        % STANDING STILL THROUGH MOST OF THE SEGMENT, IT MIGHT GET ELIMINATED
        % AS BACKGROUND. PERHAPS THERE IS A BETTER WAY TO DO THIS.

        %numbkgdframes=0; % counter

        % tot=zeros(a.Height,a.Width,3); %allocating/creating a matrix of the proper type/size
        % frameb=double(tot);

        %change in version 3.83 - using a single frame for background
        %subtraction rather than an average over the course of the video

        %disp('finding background - if frozen for >3s, hit Ctrl+C and try again')
        %fprintf('         \n')


        % using one frame every 1/6 sec is good enough for making background. 
        % start at 1 sec to avoid shaking often found at beginning of clip.

        frame=read(a,segend);
        invframe=255-frame;
        frameb=double(invframe);
        %totb=frameb;
        %numbkgdframes=numbkgdframes+1;
        %if mod(a.NumberOfFrames-k,fps)==0;
        %fprintf(1,'\b\b\b\b\b\b\b\b\b\b%10.0f',k); pause(.1)
        %end

        %fprintf('\n')
        bkgd=uint8(frameb);

        % display and save background image. This can be checked later to see
        % if any flies were eliminated.
        if verbose==1
            h=imshow(bkgd); 
            saveas(h,strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_bkgd.jpg'));
        end






        % FIND FLIES

        % takes a VideoReader object and its background, outputs a matrix containing the coordinates
        % of each fly found in each frame. X-values are in odd rows, Y-values in
        % even rows - therefore each frame takes up two rows.

        disp('finding flies')

        % Preallocate output
        flies=zeros(fix((segend-segbeg)/framerate)*2,100);

        for k=segbeg+fps+framerate:framerate:segend % for every <framerate>th 
                                                    % frame, starting after 1 sec
                                                    % to avoid jitter
            frame=read(a,k);  % get frame
            frameb=255-frame; % inverting   
            framec=frameb-bkgd; % subtracting background
            level=graythresh(framec)+.025; % changing thresholding gives different results
            framebw=im2bw(framec,level); % framebw is an array of 1's and 0's
            %framebwb=uint8(framebw); %'framebwb' is uint8 version of 'framebw' for making avi
            framebb=imopen(framebw,strel('disk',1));
            %imshow(framebb); hold on
            %ac(k).cdata=255*framebwb; % for making avi


            regions = regionprops(framebb,'basic');
                % I tried playing with the following code to make sure I was
                % detecting only flies, but 'strel' above seems to do the job.
                %    if regions(x).Area<55; %Area must be between 10 and 55
                %    if regions(x).Area>10;
                %    area = regions(x).Area

            y=1;
            for x=1:numel(regions) %  find flies, output coords into matrix 'flies'
                flies(2*((k-segbeg)/framerate)-1,y)=regions(x).Centroid(1); % putting x-value in first row
                flies(2*((k-segbeg)/framerate),y)=regions(x).Centroid(2); %  putting y-value in second row
                y=y+1;
            end


        end

        flies=round(flies); %gets rid of fractions created by "Centroid"
        flies(~flies)=nan;  %turns zeros in matrix into "Not A Number"

        %movie2avi(ac, 'ac.avi', 'compression', 'None'); % making avi

        csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_flies.csv'),flies);



        % EXTRACT GENERAL FLY DATA

        % takes array of fly positions (along horizontal axis of matrix),
        % for each frame (in pairs along vertical axis of matrix), and number of
        % vials.
        % determines number of frames in segment
        % determines height of vials by finding highest and lowest fly positions.
        % determines vial widths by finding right and leftmost flies and dividing
        % by number of vials.

        fliessize=size(flies);
        flieslength=fliessize(1); 
        numframes=fix(flieslength/2); %numframes = number of frames
        flyslots=fliessize(2); %flyslots = number of slots for flies, which is actually fixed above at 100. should maybe change this. ON SECOND THOUGHT -  I don't think it really matters what it was set at above - if it goes over, it will expand the matrix
        countedfliesinvial=zeros(numframes,numvials);
        totalheightsarray=zeros(numframes,numvials);
        numseconds=fix((numframes*framerate)/fps); %the number of seconds represented

        xvalues=zeros(numframes,flyslots);
        yvalues=zeros(numframes,flyslots);

        for k=1:numframes;    %creates a matrix of xvalues and a matrix of yvalues
            xvalues(k,:)=flies(2*k-1,:); %at this point, x and y are referring to coordinates on the image, so x is horizontal and y is vertical
            yvalues(k,:)=flies(2*k,:);
        end

        xmin=min(min(xvalues))-1;
        xmax=max(max(xvalues))+1;
        ymin=min(min(yvalues));
        ymax=max(max(yvalues));

        ymid=ymin+((ymax-ymin)/2);

        vialwidth=(xmax-xmin)/numvials;

        width(segment)=(xmax-xmin);
        heighth(segment)=(ymax-ymin);

        % OUTPUT SAMPLE FRAMES

        % saves the last image of each second of the segment, 
        % with flies and lines indicated
        % NOTE: DOES NOT SHOW FLY LOCATIONS FOR FIRST FRAME OUTPUT BECAUSE
        % FLY LOCATIONS FOR THE 1ST SECOND ARE NOT DETERMINED.
        if verbose==1
            for k=1:numseconds

            xcentroids = cat(1, flies(2*(k*fps/framerate)-1,:)); % display results on graph
            ycentroids = cat(1, flies(2*(k*fps/framerate),:));


            hold on
            h=imshow(read(a,(segbeg+k*fps)));


            line([xmin xmin],[ymin ymax])
            line([xmax xmax],[ymin ymax])
            line([xmin xmax],[ymin ymin])
            line([xmin xmax],[ymax ymax])
            for j=1:numvials
                line([xmin+j*vialwidth xmin+j*vialwidth], [ymin ymax])
            end
            line([xmin xmax],[ymid ymid])

            plot(xcentroids, ycentroids, 'r+')
            saveas(h,strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_sec',num2str(k),'.jpg'));
            if k==2
                saveas(h,strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_CHECKFRAME.jpg'));
            end

            hold off
            end

        end



        % ANALYSIS #1 - FIND AVERAGE HEIGHTS OF FLIES

        % takes fly data from above
        % returns average height of flies in vials

        for k=1:numframes %for every frame
            for j=1:flyslots % for every fly slot
                if isfinite(xvalues(k,j)) % if there's a fly in that slot
                    height=ymax-yvalues(k,j); % finds height (ymax because position is counted from top)
                    vial=ceil((xvalues(k,j)-xmin)/vialwidth); %figures out which vial the fly is in
                    totalheightsarray(k,vial)=totalheightsarray(k,vial)+height; %tallies up the total height of flies in vial in this frame
                    countedfliesinvial(k,vial)=countedfliesinvial(k,vial)+1; %countedfliesinvial is # of flies detected in vial
                    % it makes sense to do it this way, determining the number
                    % of flies in each vial anew for every frame, because not
                    % all flies will be detected in every frame. It is true
                    % that at the beginning of the segment
                    % most of the flies will be grouped at the
                    % bottom, and few flies will be detected. Doing it this way
                    % only takes into account those flies that are detected,
                    % which would tend to overestimate the average height of
                    % flies in the vial. However, this overestimate is very
                    % slight, since the values are all close to zero anyway.
                    % It is also true that at the end, many of the flies would be
                    % grouped at the top, and few of these flies would be
                    % detected, leading to underestimation of the average
                    % height. Again, this underestimate is small. 
                    % However, if the numbers of flies detected were
                    % divided by the total number of flies in the vial, the
                    % loss of flies at the top would lead to a massive
                    % underestimation of the total height. I think that to
                    % achieve the best possible dynamic range, the solution
                    % I've come up with is optimal.
                end
            end
        end


        % interesting - now I remember - I tried it with this variable and
        % eventually threw it out because I didn't like the way it worked

        countedfliesinvial(~countedfliesinvial)=1; %countedfliesinvial must be at least 1, or div by 0 error

    %    
    %     if size(varargin,2)>=1
    %         f=ones(numframes,1);
    %         inputfliesinvialmatrix=f*inputfliesinvial;
    %         fliesinvial=inputfliesinvialmatrix;
    %         avgfractionoffliesdetected=mean(countedfliesinvial./inputfliesinvialmatrix);
    %         inputfliesinvial
    %         avgfractionoffliesdetected
    %         
    %     else
    %         fliesinvial=countedfliesinvial;
    %     end


    %     fliesinvial(~fliesinvial)=1; %vialtot must be at least 1, or div by 0 error

        totalheights=(totalheightsarray)/(ymax-ymin); 
        avgheights=(totalheightsarray./countedfliesinvial)/(ymax-ymin); 
        % 'avgheights' for 'average heights of flies', 
        % which is recorded for each vial in each frame.
        % avgheights is the total heights of flies in the vial, divided by the 
        % number of flies determined to be in the vial, 
        % and expressed as a fraction of vial height.



        % CALCULATE AVERAGE HEIGHTS OVER EACH SECOND

        % avgheights is a matrix of the average height of the flies in each
        % vial for each time point evaluated. Which time points are evaluated
        % depends on the framerate. Framerate=6 indicates that every 6th frame
        % is evaluated. This routine returns a matrix with the average height of
        % the flies in each vial during the previous second. Normally,
        % there are 30 frames per second, so there would be 30/framerate frames per
        % second when passed to this routine.
        % NOTE: SINCE THE FIRST SECOND IS SKIPPED, THE AVERAGE HEIGHT RETURNED
        % FOR THE FIRST SECOND WILL ALWAYS ARTIFICIALLY BE ZERO. THIS IS NOT
        % APPARENT IN THE OUTPUT, AND COULD BE CONSIDERED A MISREPRESENTATION
        % OF THE DATA. PERHAPS THIS SHOULD BE CHANGED. I'M NOT QUITE SURE WHAT
        % TO DO WITH IT, THOUGH.
        % in versions >5, however, it doesn't really matter since the output is
        % s3-s2 (speed between seconds 2 and 3) and the first second doesn't
        % enter into it.

        onesecavg=zeros(fps/framerate,numvials); % matrix to hold the values for a single second
        onesectot=zeros(fps/framerate,numvials); % matrix to hold the values for a single second

        avgheightspersec=zeros(maxseglength,numvials);
        avgheightspersec(1,:)=NaN; % to get rid of spurious 0's

        totalheightspersec=zeros(maxseglength,numvials);
        totalheightspersec(1,:)=NaN; % to get rid of spurious 0's

        for j=2:numseconds                 % for every second, starting with the 2nd
            row=(j-1)*fps/framerate;    % row = the first entry for that second
            for k=1:fps/framerate       % k counts off the entries 
                onesecavg(k,:)=avgheights(row+k,:);
                onesectot(k,:)=totalheights(row+k,:);
            end
            avgheightspersec(j,:)=mean(onesecavg);
            totalheightspersec(j,:)=mean(onesectot);

        end

    %     csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_avgheights.csv'),avgheights)
    %     csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_avgheightspersec.csv'), avgheightspersec);
    % 
    %     csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_totalheights.csv'),avgheights)
    %     csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_seg',num2str(segment),'_totalheightspersec.csv'), avgheightspersec);
    % 


        % CONCATENATE AVGHEIGHT RESULTS ONTO PREVIOUS RESULTS MATRIX

        separator=100*segment*ones(maxseglength,1);
        separator=separator+totsecs;

        ahpsperseg=[ahpsperseg separator avgheightspersec];
        thpsperseg=[thpsperseg separator totalheightspersec];

        countedfliesinvial(1:6,:)=[]; % the first 6 frames are not evaluated, so are removed
        avgfliesinvial(segment,:)=mean(countedfliesinvial);
        stdfliesinvial(segment,:)=std(countedfliesinvial);
        maxfliesinvial(segment,:)=max(countedfliesinvial);


    end


    % OUTPUT RESULTS

    ahpsperseg(32,1)=version;
    % thpsperseg(32,1)=version; %this is being left in just in case in the future 
                                %someone wants to use the actual numbers of
                                %flies in the vials instead of the counted
                                %flies in the vials as a a divisor.


    wmedian=median(width);
    hmedian=median(heighth);


    % Eliminate segments that are out of spec. If the boxes detected are more
    % or less than 5% away from the median, they are excluded.

    excluded=0;

    for j=segments:-1:1
        if ( width(j) < ( .95 * wmedian ) || width(j) > ( 1.05 * wmedian ) || heighth(j) < ( .95 * hmedian ) || heighth(j) > ( 1.05 * hmedian ) )
            for k=1:numvials
                ahpsperseg(:,( j * (numvials+1) - k + 1 ) ) = NaN;
            end
            avgfliesinvial(j,:)=[];
            stdfliesinvial(j,:)=[];
            maxfliesinvial(j,:)=[];
            excluded=[j excluded];
        end
    end

    excludedmatrix=ones(numvials,1)*excluded;

    for j=1:numvials
        m=1;
        for k=1:segments
            if isnan(ahpsperseg(2,k*(numvials+1)))
            else    
            speeds(m,j)=ahpsperseg(3, ((k-1)*(numvials+1) + j + 1)) - ahpsperseg(2, ((k-1)*(numvials+1) + j + 1));
            m=m+1;
            end
        end
        index(1,j)=j;
    end

    avgspeeds=mean(speeds);
    stdspeeds=std(speeds);
    semspeeds=std(speeds)/sqrt(size(speeds,1)); % standard errors of the mean

    avgavgfliesinvial=mean(avgfliesinvial);
    stdavgfliesinvial=std(avgfliesinvial);
    avgstdfliesinvial=mean(stdfliesinvial);
    maxmaxfliesinvial=max(maxfliesinvial);
    avgmaxfliesinvial=mean(maxfliesinvial);
    stdmaxfliesinvial=std(maxfliesinvial);


    % other stats
    duration=ones(numvials,1)*toc;              %the first stat is the duration of the program
    segmentcount=ones(numvials,1)*segments;     %the second stat is the total number of segments


    data=[transpose([index; avgspeeds; semspeeds; stdspeeds; avgavgfliesinvial; stdavgfliesinvial; avgstdfliesinvial; maxmaxfliesinvial; avgmaxfliesinvial; stdmaxfliesinvial]) duration segmentcount excludedmatrix];
    data(1,100)=NaN;    % makes all outputs be 100 elements wide, so that they can be concatenated later

    datacell=num2cell(data);

    for j=1:numvials
        vidlabel{j,1}=dirlist(i).name;
    end

    sumcell=[vidlabel datacell];


    % csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_SUMMARY.csv'), summary);
    % 
    % csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_AVGHEIGHTS_FINAL.csv'), ahpsperseg);
    % csvwrite(strcat(dirlist(i).name,'_v',num2str(version),'_TOTALHEIGHTS_FINAL.csv'), thpsperseg);


    fprintf('\n');
    toc
    fprintf('\n');

    output=[output; sumcell];
end

headers={'VIDEO' 'VIAL' 'AVG SPEED' 'SEM SPEED' 'STD SPEED' 'AVG # FLIES' 'STD FLIES' 'AVG STD FLIES' 'MAX MAX FLIES' 'AVG MAX FLIES' 'STD MAX FLIES' 'DURATION' 'TOTAL SEGMENTS' 'EXCLUDED SEGMENTS'};
headers{1,101}=NaN;


% function cell2csv(filename,cellArray,delimiter)
% % Writes cell array content into a *.csv file.
% % 
% % CELL2CSV(filename,cellArray,delimiter)
% %
% % filename      = Name of the file to save. [ i.e. 'text.csv' ]
% % cellarray    = Name of the Cell Array where the data is in
% % delimiter = seperating sign, normally:',' (default)
% %
% % by Sylvain Fiedler, KA, 2004
% % modified by Rob Kohr, Rutgers, 2005 - changed to english and fixed delimiter
% if nargin<3
%     delimiter = ',';
% end

cellArray=[headers; output];
delimiter=',';

datei = fopen(filename,'w');
for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)

        var = eval('cellArray{z,s}');

        if size(var,1) == 0
            var = '';
        end

        if isnumeric(var) == 1
            var = num2str(var);
        end

        fprintf(datei,var);

        if s ~= size(cellArray,2)
            fprintf(datei,delimiter);
        end
    end
    fprintf(datei,'\n');
end
fclose(datei);


beep

% end cell2csv code

    
    
