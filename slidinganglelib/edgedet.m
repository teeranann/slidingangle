% White-black edge calculation.
% Teeranan Nongnual, teeranan.no@buu.ac.th
% Burapha University, THAILAND
% August 2021


%% Config parameters
marginedgediff=5;

%% Loading of data

if ispc
    if exist('pathnamm','var') && pathnamm(1)~=0
        [fnam, pathnamm, filterindex] = uigetfile('*.avi', 'Pick a video file',pathnamm);
    else
        [fnam, pathnamm, filterindex] = uigetfile('*.avi', 'Pick a video file');
    end
else
    ls *.*
    fnam=input('enter video file name ','s');
end

vidObj=VideoReader([pathnamm fnam]); %load video of tilt experiment
disp(' ')
disp('VIDEO INPUT')
disp([pathnamm fnam])
ACal= @(t) 0*t;
XCal= @(t) 0*t;
YCal= @(t) 0*t;

% read each frame in video into variable.
mov={};
k=1;
while hasFrame(vidObj)
    frame = readFrame(vidObj);
    %     mov{k}=flipud(frame(:,:,1)); %#ok<SAGROW>
    mov{k}=(frame(:,:,1));
    k = k+1;
end
camfrate=vidObj.FrameRate;
camduration=vidObj.Duration;
camframetot=round(camfrate*camduration);
camwidth=vidObj.Width;
camheight=vidObj.Height;
timest=1/camfrate*(0:camframetot-1);

fprintf('%1.0i frames\n',camframetot)
fprintf('%1.1i Hz\n',camfrate)
fprintf('%1.3f s\n',camduration)
fprintf('%1.0i x %1.0i px\n',camwidth,camheight)

plotint=round(camframetot/5);

%% Edge detection
disp(' ')
disp('EDGE CALCULATIONS')
if exist([pathnamm  fnam '-edge.mat'],'file')
    disp('Edge file detected. Loading -edge.mat file.')
    disp([pathnamm  fnam '-edge.mat'])
    load([pathnamm  fnam '-edge.mat'],'-regexp','^(?!polyf$|ellipf$|PlotPolynomialFrames$|plotpolyallframe$|polyvol$|PlotEllipticFrames$|plotellipallframe$|ellipvol$|timest$|vidObj$|zpixelcutoff$|pathnamm$|fnam$|polycal$|mov$)...')
else
    fprintf('Edge file not detected. Calculating -edge.mat from frame #1 to #%1.0i.\n',camframetot)
    % predefine matrices prior to loop
    
    EdgeCell=cell(camframetot,2); % cell structure to store detected edges
    EdgeCellLR=cell(camframetot,1);
    BaseVec=zeros(camframetot,4); % matrix to store Baseline coordinates [x0L,y0L,x0R,y0R]
    
    for n_frame=1:camframetot
        if mod(n_frame,10^(max([1 floor(log10(camframetot)-2)])))==0
            fprintf('%1.0i ', n_frame)
            if mod(n_frame,10^(max([2 floor(log10(camframetot)-2)])))==0
                fprintf('\n')
            end
        end
        if n_frame==camframetot && mod(n_frame,10^(max([1 floor(log10(camframetot)-2)])))~=0
            fprintf('\n')
        end
        im=mov{n_frame}; %extrancting image from video
        [edges, RI] = subpixelEdges(im, 25); %find edges 15

        locnotnan =~isnan(edges.x)| ~isnan(edges.y);
        edges.x=edges.x(locnotnan); edges.y=edges.y(locnotnan);
        edges.nx=edges.nx(locnotnan); edges.ny=edges.ny(locnotnan);
        edges.i0=edges.i0(locnotnan); edges.i1=edges.i1(locnotnan);
        edges.position=edges.position(locnotnan); edges.curv=edges.curv(locnotnan);
        
        % Edge for CA
        longestedge=findlongestedge(edges,size(im),marginedgediff); % select longest edge
        
        %longestedge ~= NaN
        locnotnan =~isnan(longestedge.x);
        longestedge.x=longestedge.x(locnotnan); longestedge.y=longestedge.y(locnotnan);
        longestedge.nx=longestedge.nx(locnotnan); longestedge.ny=longestedge.ny(locnotnan);
        longestedge.i0=longestedge.i0(locnotnan); longestedge.i1=longestedge.i1(locnotnan);
        longestedge.position=longestedge.position(locnotnan); longestedge.curv=longestedge.curv(locnotnan);
       
        [edgeL,edgeR]=leftrightedges(longestedge); % divide into left and right at the drop apex
        EdgeCell{n_frame,1}=edgeL; %save for later use
        EdgeCell{n_frame,2}=edgeR;
        
        fields = fieldnames(edgeL);
        edgeLR=edges;
        for i = 1:numel(fields)
            edgeLR.(fields{i})=[edgeL.(fields{i})(end:-1:1);edgeR.(fields{i})];
        end
        
        EdgeCellLR{n_frame,1}=edgeL; %save for later use
        EdgeCellLR{n_frame,2}=edgeR;
        EdgeCellLR{n_frame,1}.x=longestedge.x;
        EdgeCellLR{n_frame,1}.y=longestedge.y;
                
        EdgeCellSALR{n_frame,1}=edges;

    end
    save([pathnamm  fnam '-edge.mat'])
end
disp(' ')