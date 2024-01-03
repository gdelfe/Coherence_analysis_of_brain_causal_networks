
function smooth_arr = smooth_array(X,h) 
    
    smooth_arr = [];
    for i=1:size(X,1)
        smooth_arr = [smooth_arr; smooth(X(i,:),h)];
    end
    
end