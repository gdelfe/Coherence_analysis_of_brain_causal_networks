

function [gc_avg] = average_gc_across_monkeys(gc_mav,gc_arc,XY,YX)

XY_all = [];
YX_all = [];

for i = 1:length(gc_mav) % Maverick
    if ~isempty(gc_mav(i).(XY))
        XY_all = [XY_all; gc_mav(i).(XY)(1,2,:)];
    end
    if ~isempty(gc_mav(i).(YX))
        i
        YX_all = [YX_all; gc_mav(i).(YX)(2,1,:)];
    end
end
for i = 1:length(gc_arc) % Archie
    if ~isempty(gc_arc(i).(XY))
        XY_all = [XY_all; gc_arc(i).(XY)(1,2,:)];
    end
    if ~isempty(gc_arc(i).(YX))
        YX_all = [YX_all; gc_arc(i).(YX)(2,1,:)];
    end
end

XY_all = sq(XY_all);
YX_all = sq(YX_all);


gc_avg.(XY).mean = mean(XY_all, 1);
gc_avg.(XY).std = std(XY_all, [], 1);
n_size = size(XY_all, 1);
gc_avg.(XY).sem = gc_avg.(XY).std / sqrt(n_size);



gc_avg.(YX).mean = mean(YX_all, 1);
gc_avg.(YX).std = std(YX_all, [], 1);
n_size = size(YX_all, 1);
gc_avg.(YX).sem = gc_avg.(YX).std / sqrt(n_size);


gc_avg.freq = gc_mav.freq;


end
