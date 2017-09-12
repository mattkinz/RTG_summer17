function opt = check_align_opt(opt)

%checks the alignment options for auto align


%Written by: Toby Sanders @ Pacific Northwest National Laboratory
%Computational & Applied Mathematics Department, Univ. of South Carolina
%7/11/2014


if ~isfield(opt,'gridtype')
    opt.gridtype = 'holey';
end

if ~isfield(opt,'consistency_rating')
    opt.consistency_rating=1/5;
end

if ~isfield(opt,'window_accuracy')
    opt.window_accuracy='super';
end

if ~isfield(opt,'rescaling')
    opt.rescaling=true;
end

if opt.rescaling
    if ~isfield(opt,'rescale_value')
        opt.rescale_value=1000;
    elseif ~isnumeric(opt.rescale_value)
        opt.rescale_value=1000;
    end
end

if ~isfield(opt,'xalign_type')
    opt.xalign_type='global';
end

if ~isfield(opt,'filter')
    opt.filter=false;
    opt.filter_level=.05;
elseif ~isfield(opt,'filter_level')
    opt.filter_level=.05;
end


end
    
    