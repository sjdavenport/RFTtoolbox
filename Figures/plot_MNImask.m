MNImask = imgload('MNImask');
f = Field();
f.mask = logical(MNImask);
f.xvals = {1:size(MNImask,1),1:size(MNImask,2), 1:size(MNImask,3)};
plot_voxmf( f )