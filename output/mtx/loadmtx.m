function Data = loadmtx(name)
fp = fopen(name);
header = textscan(fgets(fp),'%s','delimiter',',');
s = sscanf(fgets(fp), '%d %d %d %d');
if (s(4) == 4)
  M = fread(fp, s(1)*s(2)*s(3), 'float');
  M = reshape(M,s(3),s(2),s(1));
else
  M = fread(fp, s(1)*s(2)*s(3), 'double');
  M = reshape(M,s(3),s(2),s(1));
Data = [header,M]; 
end