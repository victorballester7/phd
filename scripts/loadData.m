folder = "/home/victor/Desktop/PhD/src/incNSboeingGapRe1000/spod/d2_w26/";
fid = fopen(fullfile(folder, 'HistoryPoints.dat'));
points = readmatrix(fullfile(folder, 'points.txt'));
x = points(:,1);
y = points(:,2);
w = readmatrix(fullfile(folder, 'weights.txt'));

nPoints = size(points, 1);

disp("Reading data...");
data = textscan(fid, '%f %f %f %f', 'CommentStyle', '#');
disp("Data read.");
fclose(fid);

time = data{1}(1:nPoints:end);
nTime = length(time);
w = repmat(w, 2, 1)';

u = reshape(data{2}, [nPoints, nTime])';
v = reshape(data{3}, [nPoints, nTime])';


