clc;
clear all;

%%read tractography file that is in vtk polydata
Name=strcat('ThalamicRadiationR.vtk');
fid = fopen(Name,'r'); 
C= textscan(fid,'%f %f %f %f %f %f %f %f %f','HeaderLines', 5);
unordered_coords = [C{1,1},C{1,2},C{1,3},C{1,4},C{1,5},C{1,6},C{1,7},C{1,8},C{1,9}];
ID = 1;
row = 1;
Coords_ID = zeros(length(unordered_coords)*3,4);
for i=1:length(unordered_coords)
    Coords_ID(row,:) = [ID,unordered_coords(i,1:3)];
    row = row + 1;
    ID = ID + 1;
    Coords_ID(row,:) = [ID,unordered_coords(i,4:6)];
    row = row + 1;
    ID = ID + 1;
    Coords_ID(row,:) = [ID,unordered_coords(i,7:9)]; 
    row = row + 1;
    ID = ID + 1;
end
fclose(fid);
for i=1:length(Coords_ID)
    if isnan(Coords_ID(i,2))
        Coords_ID = Coords_ID(1:i-1,:);
        break
    end
end

file_path  = strcat('ThalamicRadiationR_Lines.vtk');
fileID = fopen(file_path, 'r');
lines = cell(0, 1);
while ~feof(fileID)
    current_line_data = textscan(fgetl(fileID), '%d');
    if ~isempty(current_line_data{1})
        lines{end+1, 1} = current_line_data{1};
    end
end
fclose(fileID);

Fibers=cell(0, 1);
for i=1:length(lines)
    coords_Fiber = [];
    for j=2:length(lines{i,1})
        coords_Fiber = [coords_Fiber;Coords_ID(lines{i,1}(j,1)+1,2:end)];
    end
    Fibers{end+1, 1} = coords_Fiber;
end

%% Read Surface
Name_lh_pial=strcat('rh_white.vtk');
fid = fopen(Name_lh_pial,'r'); 
C= textscan(fid,'%f %f %f %f %f %f %f %f %f','HeaderLines', 5);
Surf_lh_coords = [C{1,1},C{1,2},C{1,3},C{1,4},C{1,5},C{1,6},C{1,7},C{1,8},C{1,9}];
ID = 1;
row = 1;
Vert_ID = zeros(length(Surf_lh_coords)*3,4);
for i=1:length(Surf_lh_coords)
    Vert_ID(row,:) = [ID,Surf_lh_coords(i,1:3)];
    row = row + 1;
    ID = ID + 1;
    Vert_ID(row,:) = [ID,Surf_lh_coords(i,4:6)];
    row = row + 1;
    ID = ID + 1;
    Vert_ID(row,:) = [ID,Surf_lh_coords(i,7:9)]; 
    row = row + 1;
    ID = ID + 1;
end
fclose(fid);
for i=1:length(Vert_ID)
    if isnan(Vert_ID(i,2))
        Vert_ID = Vert_ID(1:i-1,:);
        break
    end
end

%read the faces
fid = fopen(Name_lh_pial,'r'); 
rowID = 0;
line = '';
while ~feof(fid)
    line = fgetl(fid);
    rowID = rowID + 1;
    if contains(line, 'POLYGONS')
        disp(['Row ID of the line with "POLYGONS": ', num2str(rowID)]);
        disp(['Line content: ', line]);
        break;
    end
end
fclose(fid);
fid = fopen(Name_lh_pial,'r'); 
C= textscan(fid,'%d %d %d %d','HeaderLines', rowID);
Faces = [C{1,2},C{1,3},C{1,4}]+1;
fclose(fid);

%read lookup table
fid = fopen(Name_lh_pial,'r'); 
rowID = 0;
line = '';
while ~feof(fid)
    line = fgetl(fid);
    rowID = rowID + 1;
    if contains(line, 'LOOKUP_TABLE')
        disp(['Row ID of the line with "POLYGONS": ', num2str(rowID)]);
        disp(['Line content: ', line]);
        break;
    end
end
fclose(fid);
fid = fopen(Name_lh_pial,'r'); 
C= textscan(fid,'%f %f %f %f %f %f %f %f %f','HeaderLines', rowID);
sulc = [C{1,1},C{1,2},C{1,3},C{1,4},C{1,5},C{1,6},C{1,7},C{1,8},C{1,9}];
ID = 1;
row_1 = 1;
sulc_ID = zeros(length(sulc)*3,1);
for i=1:length(sulc)
    row_2 = row_1 + 8;
    sulc_ID(row_1:row_2,:) = transpose([sulc(i,1:9)]);
    row_1 = row_2 + 1;
end
for i=1:length(sulc_ID)
    if isnan(sulc_ID(i,1))
        sulc_ID = sulc_ID(1:i-1);
        break
    end
end

%read normals
fid = fopen(Name_lh_pial,'r'); 
rowID = 0;
line = '';
while ~feof(fid)
    line = fgetl(fid);
    rowID = rowID + 1;
    if contains(line, 'NORMALS')
        disp(['Row ID of the line with "POLYGONS": ', num2str(rowID)]);
        disp(['Line content: ', line]);
        break;
    end
end
fclose(fid);
fid = fopen(Name_lh_pial,'r'); 
C= textscan(fid,'%f %f %f %f %f %f %f %f %f','HeaderLines', rowID);
Normals = [C{1,1},C{1,2},C{1,3},C{1,4},C{1,5},C{1,6},C{1,7},C{1,8},C{1,9}];
fclose(fid);

% %% find the nodes inside the surface
% Fibers_in = cell(length(Fibers),1);
% for i=1:length(Fibers)
%     inside = in_polyhedron(Faces, Vert_ID(:,2:end), Fibers{i,1});
%     Fibers_in{i,1} = Fibers{i,1}(inside,:);
%     if rem(i,100)==0
%         i
%     end
% end
%% Find the end nodes that lie on the surface
Fibers_start_end= [];
for i=1:length(Fibers)
    Fibers_start_end = [Fibers_start_end;Fibers{i,1}(1,:);Fibers{i,1}(end,:)];
end
Fibers_end_points = [];
for i=1:2:length(Fibers_start_end)
    Dist=[pdist2(Fibers_start_end(i:i+1,:),[0,0,0]),Fibers_start_end(i:i+1,:)];
    Dist = sortrows(Dist,1);
    Fibers_end_points = [Fibers_end_points;Dist(end,2:end)];
end
%% find number of fibers in sulcal and gyral
Dist = [];
sulcal_points = [];
gyral_points = [];
for i=1:length(Fibers_end_points)
    Dist=[pdist2(Vert_ID(:,2:end),Fibers_end_points(i,:)),Vert_ID(:,2:end),sulc_ID];
    Dist = sortrows(Dist,1);
    if (Dist(1,5)>=0)
        sulcal_points = [sulcal_points;Dist(1,2:4)];
    else
        gyral_points = [gyral_points;Dist(1,2:4)];
    end
end

disp(['Number of fibers on gyri: ', num2str(length(gyral_points))]);
disp(['Number of fibers on sulci: ', num2str(length(sulcal_points))]);


Name=strcat('Fibers_endpoints_R.vtk');
fileID = fopen(Name,'w');
fprintf(fileID,'# vtk DataFile Version 2.0\n');
fprintf(fileID,'Cube example\nASCII\nDATASET POLYDATA\n');
fprintf(fileID,'POINTS %d float\n',size(Fibers_end_points,1));
fprintf(fileID,'%5f %5f %5f\n',transpose(Fibers_end_points));
fclose(fileID)

%%
%areas
total_sulcal_area = 0;
total_gyri_area = 0;
total_surface_area=0;
for i=1:length(Faces)
    v1 = Vert_ID(Faces(i,1),2:end);
    v2 = Vert_ID(Faces(i,2),2:end);
    v3 = Vert_ID(Faces(i,3),2:end);
    edge1 = v2 - v1;
    edge2 = v3 - v1;
    triangle_area = 0.5 * norm(cross(edge1,edge2));
    total_surface_area = total_surface_area + triangle_area;
    
    if ((sulc_ID(Faces(i,1))>=0) && (sulc_ID(Faces(i,2))>=0) && (sulc_ID(Faces(i,3))>=0))
        total_sulcal_area=total_sulcal_area+triangle_area;
    elseif ((sulc_ID(Faces(i,1))<0) && (sulc_ID(Faces(i,2))<0) && (sulc_ID(Faces(i,3))<0))
        total_gyri_area=total_gyri_area+triangle_area;
    end
end

sulci_num_percentage = length(sulcal_points)/(length(sulcal_points) + length(gyral_points))*100;
gyri_num_percentage = length(gyral_points)/(length(sulcal_points) + length(gyral_points))*100;

FD_sulci_area = (length(sulcal_points)/total_sulcal_area);
FD_gyri_area = (length(gyral_points)/total_gyri_area);
FD_sulci_area_total = FD_sulci_area/(FD_sulci_area+FD_gyri_area)*100;
FD_gyri_area_total = FD_gyri_area/(FD_sulci_area+FD_gyri_area)*100;
