function opt = default_align_opt

%default alignment options for autoalign



%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


opt.gridtype = 'holey';
opt.consistency_rating = 1/5;
opt.window_accuracy = 'super';
opt.rescaling=true;
opt.rescale_value = 1000;
opt.xalign_type='global';
opt.filter=false;
opt.filter_level = .05;
