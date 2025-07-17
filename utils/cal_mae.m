function mae = cal_mae(aop_gt,aop)
anglediff = abs(aop-aop_gt);
minAngleDiff = min(anglediff,abs(anglediff-pi));
mae = mean(minAngleDiff(:))/pi*180;
end