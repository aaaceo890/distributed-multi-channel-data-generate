function write_info(fname,info)
% write info in txt
%   Detailed explanation goes here

fid=fopen(fname,'a');
  
fprintf(fid,'<corpus name> %s\n',info.name);

fprintf(fid,'<source info>\n');
print_coordinate(fid,info.source);

fprintf(fid,'<point-noise info>\n');
for p = 1:size(info.point,1)
    fprintf(fid,'\t<pos_%02d>',p);
    print_coordinate(fid,info.point(p,:));
end

fprintf(fid,'<mic info>\n');
for m = 1:size(info.mic,1)
    fprintf(fid,'\t<mic_%02d>',m);
    print_coordinate(fid,info.mic(m,:));
    fprintf(fid,[repmat('\t',1,4),'<snr> %.2f\n'],info.snr(m));
end
  
fprintf(fid,'<T60>\n');
for m = 1:size(info.mic,1)
    fprintf(fid,'\t<mic_%02d>',m);
    fprintf(fid,'\t<T60> %.4f\n',info.T60(m));        
end 


fprintf(fid,[repmat('=',1,60),'\n']);

fclose(fid);
end

function print_coordinate(fid,x)
fprintf(fid,'\t<x> %.2f <y> %.2f <z>%.2f\n',x(1),x(2),x(3));
end
