% generate a set of transforms that can be used to test a series of
% transforms for two set of objects

NUM_ELEMS = 8192;
FILE_NAME = 'tasks0.csv';
Mu = [500, 500, 50];
Range = [2000, 2000, 50];

fileID = fopen(FILE_NAME, 'w');
fprintf(fileID, 'i1 i2 dist\n');

for i = 1:NUM_ELEMS
    val = (rand(1,3) - 0.5)/0.5.*Range + Mu;
    fprintf(fileID, '%d ', int32(val(1)));
    fprintf(fileID, '%d ', int32(val(2)));
    fprintf(fileID, '%.7f\n', val(3));
end

fclose(fileID);
