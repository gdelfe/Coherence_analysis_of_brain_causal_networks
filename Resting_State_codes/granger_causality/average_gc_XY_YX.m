% Average Granger results (along frequency) across sessions
% Return a structure with mean, std, and sem for both GC XY and YX

function [gc_stat] = average_gc_XY_YX(gc,XY,YX)
    XY_all = [];
    YX_all = [];
    
    for i = 1:length(gc)
        if ~isempty(gc(i).(XY))
            XY_all = [XY_all; gc(i).(XY)(1,2,:)]; 
        end
        if ~isempty(gc(i).(YX))
            YX_all = [YX_all; gc(i).(YX)(2,1,:)]; 
        end
    end
    XY_all = sq(XY_all);
    YX_all = sq(YX_all);

    if ~isempty(XY_all)
        gc_stat.(XY).mean = mean(XY_all, 1);
        gc_stat.(XY).std = std(XY_all, [], 1);
        n_size = size(XY_all, 1);
        gc_stat.(XY).sem = gc_stat.(XY).std / sqrt(n_size);
    else
        gc_stat.(XY).mean = [];
        gc_stat.(XY).std = [];
        gc_stat.(XY).sem = [];
    end
    
    if ~isempty(YX_all)
        gc_stat.(YX).mean = mean(YX_all, 1);
        gc_stat.(YX).std = std(YX_all, [], 1);
        n_size = size(YX_all, 1);
        gc_stat.(YX).sem = gc_stat.(YX).std / sqrt(n_size);
    else
        gc_stat.(YX).mean = [];
        gc_stat.(YX).std = [];
        gc_stat.(YX).sem = [];
    end

    gc_stat.freq = gc.freq;
    

end
