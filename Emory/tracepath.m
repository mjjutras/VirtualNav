clear all
for sessions = 1107:1107

    
    file = ['C:\Python27\data\Giz\session_' num2str(sessions) '\log.txt'];
    fid = fopen(file);
    line_count = 1;
    l = 1;
    while l>0
       l = fgets(fid);
       if l>0
           txt{line_count,1} = l;
           line_count = line_count+1;
       end
    end

    a = 1;
    for t = 1:length(txt)
        if length(txt{t}) > 52
        if strcmp(txt{t}(46:53), 'LPoint3f')
            testa{a, 1} = txt{t}(54:end);
            a = a+1;
        end
        end
    end

    for i = 1:length(testa)
        b = find(testa{i}==',');
        xy(i, 1) = str2num(testa{i}(2:b(1)-1));
        xy(i, 2) = str2num(testa{i}(b(1)+2:b(2)-1));
    end


    sessions
    figure;
    xlim([-12 12]); ylim([-16 16]);
    hold all
    comet(xy(:, 1), xy(:,2)); figure(gcf);
    clear xy txt testa
end