function longestedge=findlongestedge(edges,imsize,margin)
% This function takes all edges discovered by the partial area edge
% detection algorithm, joins edges that are closer to eachother than the
% margin. It returns the longest edge. Input is
% edges in format
%       edges.x         xcoordinates, subpixel
%       edges.y         ycoordinates, subpixel
%       edges.porison   linear index of x,y position in pixel
% imsize in format
%       [row,col]
% margin, pixel distance between edges that should be considered part of
%                                                       same contour.


% First we determine the amount of seperate edges found spaced by more than
% 1 pixel. This is done in pixel resolution with bwlabem matlab function


logical=zeros(imsize);
logical(edges.position)=1; % create binary image with wite pixels at edges and black where there is no edges.
[L,NumberOfboundaries]=bwlabel(logical,8); % find and label seperate edges

%Load the seperate boundaries into a cell structure that have dimension
%{number of boundaries, 1}

boundary=cell(NumberOfboundaries,1); 
for ii=1:NumberOfboundaries
    [row,col]=find(L==ii);
    boundary{ii}=[row,col];
end


% Calculate the minimum distance between all points in one boudnary and all
% other points in another boundary. This is to be used for joining
% boundaries that are closer than the threshold set in margin input
% variable. The following loop creates a matrix similar to the output from
% pdist2 but wuith the minimum distance between points boundary ii and ss
% as the value in distmatrix(ii,ss) and distmatrix(ss,ii)

distmatrix=NaN(NumberOfboundaries);

for ii=1:NumberOfboundaries
    if ii+1<=NumberOfboundaries
    for ss=ii+1:NumberOfboundaries
    dist = pdist2(boundary{ii},boundary{ss});
    distmatrix(ii,ss)=min(min(dist));
    distmatrix(ss,ii)=distmatrix(ii,ss);
    end 
    end
end

% Now we have calculated all relevant distances between boundaries found and
% can start joining boudarys that are closely spaced. We start by giving
% each boundary a seperate index and the giving changing the index of
% closely space boundaries to the index of their neighbour

BoundaryIndex=1:NumberOfboundaries;
for ii=1:NumberOfboundaries
    for ss=1:NumberOfboundaries
        if distmatrix(ii,ss)<margin
            BoundaryIndex(BoundaryIndex==BoundaryIndex(ss))=BoundaryIndex(ii);
        end
   end
end

%Now there is several boundaries with same index and we combine the
%boundaries with same index into same cell

NewIndex=unique(BoundaryIndex); % the remaining boundary indexes after renaming
NumberOfNewIndexes=length(NewIndex); %number of remaining boundaries
CombinedBoundary=cell(1,NumberOfNewIndexes);


for ii=1:NumberOfboundaries
    index=find(NewIndex==BoundaryIndex(ii));
    CombinedBoundary{index}=[CombinedBoundary{index};boundary{ii}]; %joining boundaries of same index
end


% Calculating length the joined boundaries to find the longest one 
 boundarylength=zeros(1,NumberOfNewIndexes);
for ii=1:NumberOfNewIndexes
     boundarylength(ii)=length(CombinedBoundary{ii});
end

[~,index]=max(boundarylength);
mainedge=CombinedBoundary{index}; %select the longest of boundaries

[I,J] = ind2sub(imsize,edges.position); %calculate the in

pts=[I,J]; 
[~,loc] = ismember(mainedge,pts,'rows'); %calculate the location of the longest edge in the original edge structure

fields = fieldnames(edges);
for i = 1:numel(fields)
longestedge.(fields{i})=edges.(fields{i})(loc); %asign all properties of the input edge structure to the new longestedge structure
end


