function visualize(file1, file2, id1, id2)
f1 = fopen(file1, 'r');
f2 = fopen(file2, 'r');
header1 = textscan(f1, '%s%f%s%f', 1, 'delimiter', ' ');
header2 = textscan(f2, '%s%f%s%f', 1, 'delimiter', ' ');

bv_fmt =    '';
tri_fmt = '';
for i = 1:19
    bv_fmt = strcat(bv_fmt, '%f');
end
for i = 1:9
    tri_fmt = strcat(tri_fmt, '%f');
end

bv1 = cell2mat(textscan(f1, bv_fmt , header1{2}, 'delimiter', ' '));
bv2 = cell2mat(textscan(f2, bv_fmt , header2{2}, 'delimiter', ' '));

tri1 = cell2mat(textscan(f1, tri_fmt , header1{4}, 'delimiter', ' '));
tri2 = cell2mat(textscan(f2, tri_fmt , header2{4}, 'delimiter', ' '));

close all;
hold on;
plotBV(bv1(id1 + 1, :), 'r', tri1);
plotBV(bv2(id2 + 1, :), 'b', tri2);
hold off;
axis equal;

function plotBV(bv, col, tri)
p = zeros(5,3);
p(1,:) = bv(10:12);
l1 = bv(13);
l2 = bv(14);
p(2,:) = p(1,:) + l1*bv(1:3);
p(4,:) = p(1,:) + l2*bv(4:6);
p(3,:) = p(2,:) + p(4,:) - p(1,:);
p(5,:) = p(1,:);
plot(p(:,1), p(:,2), strcat(col,'-*'))

if(bv(19) >= 0)
    plotTri(tri(bv(19) + 1, :), col);
end

function plotTri(tri, col)
t = zeros(4,3);
t(1:3,:) = reshape(tri, [3,3])';
t(4,:) = t(1,:);
plot(t(:,1), t(:,2), strcat(col,'--+'))
