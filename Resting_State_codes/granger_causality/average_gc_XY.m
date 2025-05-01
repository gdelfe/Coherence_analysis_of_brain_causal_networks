
function gc_avg = average_gc_XY(gc_mav,gc_arc,XY_all, YX_all)

XY_all = [gc_mav.(XY_all); gc_arc.(XY_all)];
YX_all = [gc_mav.(YX_all); gc_arc.(YX_all)];


gc_avg.XY_mean = mean(XY_all, 1);
gc_avg.XY_std = std(XY_all, [], 1);
n_size = size(XY_all, 1);
gc_avg.XY_sem = gc_avg.XY_std / sqrt(n_size);



gc_avg.YX_mean = mean(YX_all, 1);
gc_avg.YX_std = std(YX_all, [], 1);
n_size = size(YX_all, 1);
gc_avg.YX_sem = gc_avg.YX_std / sqrt(n_size);


gc_avg.freq = gc_mav.freq;

end
