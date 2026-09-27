clear;
clc;
%start timer (just a check dw)
tic
%opening all .dcmfiles

%%copy the path into here to run the program
folderPath = ('C:\Users\ameli\OneDrive\Desktop\dcmfiles\c2');
filePattern = fullfile(folderPath, '*.dcm');
dcmFiles = dir(filePattern);

numFiles = length(dcmFiles);

fprintf('Found %d DICOM files.\n', numFiles);

allImages = cell(1, numFiles);
allInfo = cell(1, numFiles);

for i = 1:numFiles
    
    fullFileName = fullfile(dcmFiles(i).folder, dcmFiles(i).name);
    
    fprintf('Reading file %d of %d: %s\n', ...
        i, numFiles, dcmFiles(i).name);
    
    allImages{i} = dicomread(fullFileName);
    allInfo{i} = dicominfo(fullFileName);
    
end

numSlices = numFiles;

%convert all images to houndsville (dawgsville)

for i = 1:numFiles
    
    pixelData = double(allImages{i});
    
    slope = allInfo{i}.RescaleSlope;
    intercept = allInfo{i}.RescaleIntercept;
    
    allImages{i} = pixelData * slope + intercept;
    
end

%build volume (x,y, slice)

[rows, columns] = size(allImages{1});

volume = zeros(rows, columns, numSlices);

for n = 1:numFiles
    volume(:,:,n) = allImages{n};
end

%see if sclices are in order 

zPositions = zeros(numFiles,1);

for i = 1:numFiles
    
    position = allInfo{i}.ImagePositionPatient;
    zPositions(i) = position(3);
    
end

%display 3D image (.stl)

figure;
plot(1:numFiles, zPositions, 'o-');

xlabel('File Number');
ylabel('Z Position (mm)');
title('DICOM Slice Positions');

grid on;

%shows the image stack

figure;
sliceViewer(volume);

%%smoothing segments out (CHANGE THIS)

yes = volume > 1;

% Remove small objects
yes = bwareaopen(yes,100);

% Fill holes
yes = imfill(yes,'holes');

%fill small gaps

se = strel('sphere',2);
yes = imclose(yes,se);

%fill holes again after closing gaps

yes = imfill(yes,'holes');

%make da 3d surface

[F,V] = isosurface(yes,0.5);

%show 3d da surface

figure;

patch('Faces',F,'Vertices',V, ...
      'FaceColor',[0.8 0.8 0.8], ...
      'EdgeColor','none');

axis equal;

xlabel('X');
ylabel('Y');
zlabel('Z');

title('3D Bone Surface');

camlight;
lighting gouraud;

%export as .stl

TR = triangulation(F,V);

outputFile = fullfile(folderPath,'sacrum.stl');
%export .stl for before smoothing
stlwrite(TR,outputFile);

fprintf('STL file created:\n%s\n',outputFile);


%mesh smoothing (CHANGE THIS for more smooth)

numIterations = 10;
lambda = 0.2;  %.2 is sweet spot

numVertices = size(V,1);

%create adjacency matrix

A = sparse(numVertices,numVertices);

for i = 1:size(F,1)
    
    v1 = F(i,1);
    v2 = F(i,2);
    v3 = F(i,3);
    
    A(v1,v2) = 1;
    A(v2,v1) = 1;
    
    A(v1,v3) = 1;
    A(v3,v1) = 1;
    
    A(v2,v3) = 1;
    A(v3,v2) = 1;
    
end

%find number of neighbors

degree = sum(A,2);

%smooth the mesh

V_smooth = V;

for iteration = 1:numIterations
    
    %find average position of neighboring vertices
    
    neighborAverage = A * V_smooth;
    
    validVertices = degree > 0;
    
    neighborAverage(validVertices,:) = ...
        neighborAverage(validVertices,:) ./ ...
        degree(validVertices);
    
    %move vertices toward neighboring vertices
    
    V_smooth(validVertices,:) = ...
        (1-lambda) * V_smooth(validVertices,:) + ...
        lambda * neighborAverage(validVertices,:);
    
end

%show smoothed surface

figure;

patch('Faces',F,'Vertices',V_smooth, ...
      'FaceColor',[0.8 0.8 0.8], ...
      'EdgeColor','none');

axis equal;

xlabel('X');
ylabel('Y');
zlabel('Z');

title('Mesh Smoothed Sacrum');

camlight;
lighting gouraud;

%export smoothed STL

TR_smooth = triangulation(F,V_smooth);

smoothOutputFile = fullfile(folderPath, ...
    'c2_mesh_smoothed.stl');
%export smoothed STL
stlwrite(TR_smooth,smoothOutputFile);

fprintf('Smoothed STL file created:\n%s\n',smoothOutputFile);


%comparison graph

figure;

%original

subplot(1,2,1);

patch('Faces',F,'Vertices',V, ...
      'FaceColor',[0.8 0.8 0.8], ...
      'EdgeColor','none');

axis equal;

xlabel('X');
ylabel('Y');
zlabel('Z');

title('Original');

camlight;
lighting gouraud;


%smoothed

subplot(1,2,2);

patch('Faces',F,'Vertices',V_smooth, ...
      'FaceColor',[0.8 0.8 0.8], ...
      'EdgeColor','none');

axis equal;

xlabel('X');
ylabel('Y');
zlabel('Z');

title('Smoothed');

camlight;
lighting gouraud;
toc
