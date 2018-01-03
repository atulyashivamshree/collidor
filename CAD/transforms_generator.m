% generate a set of transforms that can be used to test a series of
% transforms for two set of objects

NUM_ELEMS = 10;
FILE_NAME = 'transforms_set0.csv';
Mu = [0, 0, 0, 0, 0, 0];
Range = [pi, pi, pi, 100, 100, 100]* 0 ;

fileID = fopen(FILE_NAME, 'w');
fprintf(fileID, 'Roll Pitch Yaw X Y Z\n');

for i = 1:NUM_ELEMS
    val = (rand(1,6) - 0.5)/0.5.*Range + Mu;
    fprintf(fileID, '%.7f %.7f %.7f %.7f %.7f %.7f\n', val);
end

fclose(fileID);
