%% Parameters:
%%  vr  - initialized control structure
%%  res - results cell vector
%%  per - parameters setup
%%  opt - cell vector of options, following options are available:
%%          {'mode',hist_n} - mode selector (first option, required):
%%            'norm' - dependence along X-vector
%%            'mean' - dependence along X-vector with means along Y-vector
%%            'hist' - histogram along X-vector with hist_n bins
%%            hist_n - use only for histogram mode
%%          {'var_name'} - result struct element on graph Y-axis (second option, required)
%%          {'logx'} - enable x-log scale
%%          {'logy'} - enable y-log scale
%%          {'stdbar',exp_k} - plot error bars, applies only for 'mean' mode, 
%%            the bars are generated: exp_k*std(Y-vector) 
%%          {'lim','la_var_name','lb_var_name'} - plot absolute limits, from variables 'la_var_name' and 'lb_var_name'
%%          {'lim','l_var_name',exp_k} - plot relative limits, the limit is calculated from 'l_var_name',
%%            for 'diff' option: limit = +/- exp_k*l_var_name
%%            otherwise:         limit = var_name +/- exp_k*l_var_name
%%          {'diff','d_var_name'} - difference mode, plotted values: var_name - d_var_name
%%          {'ref','r_var_name'} - plots reference value from 'r_var_name'
%%  varargin - setup to querry data from results vector:
%%           - if nothing defined X-vector is forst vector paramter in 'par' struct
%%           - first paramter is X-vector name
%%           - second paramter for 'mean' mode is Y-vector name
%%           - following parameters indexes vector elements of par struct, expected format:
%%             'par_name_1',index_1,'par_name_2',index_2,...
%%
%% Example data:
%%  Simulator setup:
%%   par.A = 1;
%%   par.B = 1:1000;
%%   par.C = [0.1 0.2 0.5 1]; 
%%   par.D = [3 4 5];
%%   par.E = 0;
%%
%%  Simulation results struct (all scalar elements):
%%   res.k2
%%   res.uk2
%%   res.k2_real
%%
%% Example function paramters:
%%  Plot simple results dependence along B vector, for par.C(1), par.D(1), with y-logscale, 
%%  use 'res{}.k2' for y-asix:
%%   ...,{'norm','k2','logy'}
%%   ...,{'norm','k2','logy'},'B'
%%
%%  Plot simple results dependence along B vector, for par.C(3), par.D(2), with y-logscale,
%%  use 'res{}.k2' for y-axis:
%%   ...,{'norm','k2','logy'},'B','C',3,'D',2
%%
%%  The same with limits from 'res.uk2' variable expanded by k = 2:
%%   ...,{'norm','k2','logy','lim','uk2',2},'B','C',3,'D',2
%%
%%  Plot results dependence along par.C with mean along par.B, for par.D(3), 
%%  use 'res{}.k2' for y-axis:
%%   ...,{'mean','k2'},'C','B','D',3
%%
%%  The same but plot difference: 'res{}.k2 - res{}.k2_real':
%%   ...,{'mean','k2','diff','k2_real'},'C','B','D',3
%%
%%  The same, include 3*std(res{}.k2) along B vector
%%   ...,{'mean','k2','diff','k2_real','stdbar',3},'C','B','D',3
%%
%%  Plot histogram with 20 bins along B vector for par.C(1) and par.D(2), 
%%  use res{}.k2 for x-axis:
%%   ...,{'hist',20,'k2'},'B','D',2

function [] = var_plot_results_vect_ieee(vr,res,par,opt,varargin)
    
  if(nargin()<4)
    % not enough paramters
    error(['Not enough parameters!']);
  endif
  
  if(~iscell(opt))
    % mode is not cell
    error(['Not enough parameters!']);
  endif
  
  % defult setup
  mode = 'norm';
  hist_n = 0;
  logx = 0;
  logy = 0;
  std_y = 0;
  lim_y = 0;
  lim_y_name = '';
  lim_y_name_b = '';
  diff_v_name = '';
  ref_v_name = '';
    
  %% parse options
  if(numel(opt))
    % identify mode
    modes = {'norm','mean','hist'};
    mid = find(strcmp(opt{1},modes),1,'first');
    if(~numel(mid))
      % undefined mode
      error(['Mode ''' opt{1} ''' undefined!']);
    endif
    mode = modes{mid};
    opt = {opt{2:end}};
    
    % get paramter for 'hist' mode
    if(strcmp(mode,'hist') && (numel(opt)<1 || ~isscalar(opt{1})))
      error('Missing paramter for ''hist'' option!');
    elseif(strcmp(mode,'hist'))
      hist_n = opt{1};
      opt = {opt{2:end}};
    endif
        
    if(numel(opt)>=1)
      % identify y-axis parameter name 
      par_name = opt{1};
      if(~isfield(res{1},par_name))
        error(['Element ''' par_name ''' is not member of results structure!']);
      endif
    
      if(numel(opt)>=1)    
        % look for 'logx' option
        oid = find(strcmp(opt,'logx'),1,'first');
        if(numel(oid))
          opt = {opt{1:oid-1} opt{oid+1:end}};
          logx = 1;  
        endif
            
        if(numel(opt)>=1)
          % look for 'logy' option
          oid = find(strcmp(opt,'logy'),1,'first');
          if(numel(oid))
            opt = {opt{1:oid-1} opt{oid+1:end}};
            logy = 1;  
          endif
        
          if(numel(opt)>=2)
            % look for 'stdbar' option {'stdbar', expansion_coef}
            oid = find(strcmp(opt,'stdbar'),1,'first');
            if(numel(oid))
              if(oid==numel(opt) || ~isscalar(opt{oid+1}))
                error('Invalid or missing value for ''stdbar'' option!');
              endif
              std_y = opt{oid+1};
              opt = {opt{1:oid-1} opt{oid+2:end}};          
            endif
            
            if(numel(opt)>=3)
              % look for 'lim' option {'lim', 'var_name', expansion_coef}
              % or for asymetrical limits {'lim', 'var_a_name', 'var_b_name'}
              oid = find(strcmp(opt,'lim'),1,'first');
              if(numel(oid))
                if(numel(opt)-oid+1<3 || ~ischar(opt{oid+1}) || (~isscalar(opt{oid+2}) && ~ischar(opt{oid+2})))
                  error('Invalid or missing parameters for ''lim'' option!');
                endif
                if(ischar(opt{oid+2}))
                  lim_y = 1;
                  lim_y_name_b = opt{oid+2};
                else
                  lim_y = opt{oid+2};
                endif              
                lim_y_name = opt{oid+1};
                opt = {opt{1:oid-1} opt{oid+3:end}};
              endif
              
              if(numel(opt)>=2)
                % look for 'diff' option {'diff','var_name'}
                oid = find(strcmp(opt,'diff'),1,'first');
                if(numel(oid))
                  if(numel(opt)-oid+1<2 || ~ischar(opt{oid+1}))
                    error('Invalid or missing parameters for ''diff'' option!');
                  endif
                  diff_v_name = opt{oid+1};
                  opt = {opt{1:oid-1} opt{oid+2:end}};              
                endif
                
                if(numel(opt)>=2)
                  % look for 'ref' option {'ref','var_name'}
                  oid = find(strcmp(opt,'ref'),1,'first');
                  if(numel(oid))
                    if(numel(opt)-oid+1<2 || ~ischar(opt{oid+1}))
                      error('Invalid or missing parameters for ''ref'' option!');
                    endif
                    ref_v_name = opt{oid+1};
                    opt = {opt{1:oid-1} opt{oid+2:end}};            
                  endif
                endif
              endif      
            endif               
          endif
        endif
      endif    
    endif
  endif
             
  if(strcmp(mode,'mean') && length(varargin)>=1)
    %%%% dependence with mean+std %%%%
  
    %% get result matrix
    [v,x_val,y_val,x_name,y_name,info_str] = var_get_results_mat(res,par,vr,varargin{:});
    
    %% add mean-axis info to info_str
    mn_str = ['mean(' y_name ')'];
    if(numel(info_str))
      info_str = [mn_str ', ' info_str];
    else
      info_str = mn_str;
    endif 
            
    %% get main element
    % found?
    if(~isfield(v,par_name))
      error(['Element ''' par_name ''' is not member of results structure!']);
    endif
    % get its matrix    
    r_var_v = getfield(v,par_name);
    % calculate mean+std
    r_var = mean(r_var_v);
    r_var_std = std(r_var_v)*std_y;
    
    %% get 'lim' element
    % found?
    if(lim_y && ~isfield(v,lim_y_name))
      error(['Element ''' lim_y_name ''' is not member of results structure!']);
    elseif(lim_y && numel(lim_y_name_b) && ~isfield(v,lim_y_name_b))
      error(['Element ''' lim_y_name_b ''' is not member of results structure!']);
    elseif(lim_y)
      
      % calculate mean
      if(numel(lim_y_name_b))
        % asymetric limits
        l_var_v = getfield(v,lim_y_name);
        l_var_a = mean(l_var_v) - r_var;
        l_var_v = getfield(v,lim_y_name_b);
        l_var_b = mean(l_var_v) - r_var;
      else
        % symetric limits
        l_var_v = getfield(v,lim_y_name);
        l_var_a = -mean(l_var_v)*lim_y;
        l_var_b = +mean(l_var_v)*lim_y;
      endif
      
        
    endif
    
    %% get 'diff' element
    % found?
    if(numel(diff_v_name) && ~isfield(v,diff_v_name))
      error(['Element ''' diff_v_name ''' is not member of results structure!']);
    elseif(numel(diff_v_name))
      % yaha, get its matrix
      d_var_v = getfield(v,diff_v_name);
      % calculate mean
      d_var = mean(d_var_v);
    endif
    
    %% get 'ref' element
    % found?
    if(numel(ref_v_name) && ~isfield(v,ref_v_name))
      error(['Element ''' diff_v_name ''' is not member of results structure!']);
    elseif(numel(ref_v_name))
      % yaha, get its vector
      ref_var = mean(getfield(v,ref_v_name));
    endif
      
  else
    %%%% single value dependence %%%%
    
    %% get vector of results
    [r,v,x_val,x_name,info_str] = var_get_results_vect(res,par,vr,varargin{:});
    
    %% get main element
    % found?
    if(~isfield(v,par_name))
      error(['Element ''' par_name ''' is not member of results structure!']);
    endif
    % get its vector    
    r_var = getfield(v,par_name);
    
    %% get 'lim' element
    % found?
    if(lim_y && ~isfield(v,lim_y_name))
      error(['Element ''' lim_y_name ''' is not member of results structure!']);
    elseif(lim_y && numel(lim_y_name_b) && ~isfield(v,lim_y_name_b))
      error(['Element ''' lim_y_name_b ''' is not member of results structure!']);
    elseif(lim_y)
      % get its vector
      l_var = getfield(v,lim_y_name)*lim_y;
      
      % calculate mean
      if(numel(lim_y_name_b))
        % asymetric limits
        l_var_a = getfield(v,lim_y_name) - r_var;
        l_var_b = getfield(v,lim_y_name_b) - r_var;
      else
        % symetric limits
        l_var_v = getfield(v,lim_y_name);
        l_var_b = getfield(v,lim_y_name)*lim_y;
        l_var_a = -l_var_b;        
      endif
      
    endif
        
    %% get 'diff' element
    % found?
    if(numel(diff_v_name) && ~isfield(v,diff_v_name))
      error(['Element ''' diff_v_name ''' is not member of results structure!']);
    elseif(numel(diff_v_name))
      % yaha, get its vector
      d_var = getfield(v,diff_v_name);
    endif
    
    %% get 'ref' element
    % found?
    if(numel(ref_v_name) && ~isfield(v,ref_v_name))
      error(['Element ''' diff_v_name ''' is not member of results structure!']);
    elseif(numel(ref_v_name))
      % yaha, get its vector
      ref_var = getfield(v,ref_v_name);
    endif
        
  endif
  
   
  %% create Y plot vector
  if(numel(diff_v_name))
    % difference plot mode
    y = r_var - d_var;
  else
    % direct plot mode
    y = r_var;  
  endif
  
   
  figure;
  if(strcmp(mode,'hist'))
    %%%% plot histogram %%%%
    
    %% plot histogram
    hist(y,hist_n);
    
    %% add some formatting
    hold off;  
    box on;
    grid on;
    if(numel(info_str))
      info_str = [', ' info_str];
    endif
    tit = ['Histogram for variable ''' x_name '''' info_str];
    %title(str_insert_escapes(tit),'FontWeight','bold');
    ylabel('N [-]','FontWeight','bold');
    if(numel(diff_v_name))
      % 'diff' mode
      xlabel(['\Delta(' str_insert_escapes(par_name) ')']);      
    else
      % direct plot mode
      xlabel([str_insert_escapes(par_name)]);
    endif
  
  else
    
    %%%% plot graph %%%%   
    %% plot main dependence
    if(logx)
      if(logy)
        loglog(x_val,y,'+--','LineWidth',3);
      else
        semilogx(x_val,y,'+--','LineWidth',3);
      endif
    else
      if(logy)
        semilogy(x_val,y,'+--','LineWidth',3);
      else
        plot(x_val,y,'+--','LineWidth',3);
      endif
    endif
    hold on;
    
    %% plot limits 
    if(lim_y)
      ofs = y;
      if(numel(diff_v_name))
        ofs = 0;
      endif
      % plot upper and lower limit
      plot(x_val,ofs+l_var_a,'k','LineWidth',1.5);
      plot(x_val,ofs+l_var_b,'k','LineWidth',1.5);
    endif
    
    %% plot reference
    if(numel(ref_v_name))
      plot(x_val,ref_var,'color',[0.8 0.0 0.0],'LineWidth',4);
    endif
    
    %% plot failed test markers (only for diff mode)
    if(numel(diff_v_name) && lim_y)
      eid = find(y>l_var_b | y<l_var_a);
      plot(x_val(eid),y(eid),'ro','MarkerSize',4,'LineWidth',2);
    endif
    
    %% plot error bars (only for 'mean' mode)
    if(std_y)
      errb = [y + r_var_std;y - r_var_std];
      %plot([x_val;x_val],errb,'b-','LineWidth',2);
      h = errorbar(x_val,y,r_var_std,'~.');
      set(h,'Marker','none');  
      set(h,'LineWidth',2);
    endif
    
    %% some formatting
    hold off;  
    box on;
    grid on;
    if(numel(info_str))
      info_str = [', ' info_str];
    endif
    tit = ['Simulator output for variable ''' x_name '''' info_str];
    %title(str_insert_escapes(tit),'FontWeight','bold');
    xlabel(str_insert_escapes(x_name));
    if(numel(diff_v_name))
      % 'diff' mode
      ylabel(['\Delta(' str_insert_escapes(par_name) ')']);
      legend(['\Delta(' str_insert_escapes(par_name) ')'],['lim(\Delta(' str_insert_escapes(par_name) '))']);
    else
      % direct plot mode
      ylabel([str_insert_escapes(par_name)]);
      legend([str_insert_escapes(par_name)],['lim(' str_insert_escapes(par_name) ')']);   
    endif
  endif

endfunction
