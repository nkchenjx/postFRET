function save4vbFRET(data, filename)
a = data(2:3,:);
a = a';
dlmwrite(filename, a, '\t')
end