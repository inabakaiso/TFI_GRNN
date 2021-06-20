function Data = Decoder(y, g, mem)
%*****************************************************************
% Decoder is convolution decoder using soft dicision Viterbi algorithm.
%   Data = Decoder(Y, G) decodes the received code Y with the transfer
%   function matrix G. The code Y should be binary.
%
%   Data = Decoder(Y, G, MEM) keeps the trace memory length no greater
%   than MEM.
%
%   INPUTS
%   =======
%   y = Coded data
%   g = Generate vector
%   mem=Trace memory
%
%   OUTPUT
%   ========
%   Data= Decoded data
%
%  Copyright (c) Chang-Jun Ahn 2000 (jun@sasase.ics.keio.ac.jp) 
%******************************************************************

if nargin < 2
    error('Not enough input variables for vit_h.');
elseif nargin < 3
    mem = 0;
end;

[n, k] = size(g);       % coding rate=1/n, constraint length=k
m = k-1;                % number of shift registor
n_shift_state = 2^m;    % number of shift registor states in trellis
n_out_state = 2^n;      % number of encoder output states

if min(size(y)) ~= 1
     error('Input must be a vector.');
end;

[y, nz] = s2p(y(:), n); % serial to n parallel
if nz ~= 0              % check the code length
     error('Code length does not match the coding rate.');
end;

len_total = size(y,1);  % length of information + termination bits
len_x = len_total - m;  % length of information bits
if len_x <= 0
     error('Code length is not enough to decode.');
end;

%****************************************************************
% The case of unlimited memory length.
%****************************************************************
if (mem == 0) | (len_total <= mem)
     mem = len_total;
end;

for i=1:n_out_state;
     out_pat(i,:) = d2b(i-1,n);      % encoder output pattern [0 1]
end;
out_pat=(-1).^out_pat;              % encoder output pattern [-1, 1]

for i = 1:n_shift_state
     shift_state(i,:) = d2b(i-1,m);  % converts state number to binary
     [out_0, state_0] = encode_bit(g, 0, shift_state(i,:));
     [out_1, state_1] = encode_bit(g, 1, shift_state(i,:));
     % encoder output patter number [(input 0) (input 1)]
     out_pat_n(i,:) = [(b2d(out_0) + 1) (b2d(out_1) + 1)];
     % transition matrix [(input 0) (input 1)]
     transition(i,:) = [(b2d(state_0)+1) (b2d(state_1)+1)];
end;
clear out_0 out_1 state_0 state_1 shift_state;

%****************************************************************
% initialize surviving path memory and path metric matrices
%****************************************************************
path_mem = NaN*ones(n_shift_state,mem);
path_metric = NaN*ones(n_shift_state,1);

%****************************************************************
% determine branch metric and path memory matrices at time 1
%****************************************************************
for i=0:1           % Input is 0 or 1
     next_state = transition(1,i+1);
     branch_metric = sum((out_pat(out_pat_n(1,i+1),:) - y(1,:)).^2);
     path_metric(next_state) = branch_metric;
     path_mem(next_state,1) = i;
end;

path_mem_tmp = path_mem;

if mem == len_total % The case of unlimited memory length
     for time = 2:len_total
         path_met_tmp = inf*ones(n_shift_state,1);
         % calculate branch metrics against every output pattern
         for i = 1:size(out_pat,1)
             branch_metric(i) = sum((out_pat(i,:) - y(time,:)).^2);
         end;
         s=1;
         while s <= n_shift_state
             if path_metric(s) == NaN
                 s = s+1;
             else
                 for i = 0:1     % Input is 0 or 1
                     next_state = transition(s,i+1);
                     branch(next_state) = branch_metric(out_pat_n(s,i+1));
                     if path_met_tmp(next_state) >=path_metric(s)+branch(next_state);
                         path_met_tmp(next_state) =path_metric(s)+branch(next_state);
                         path_mem_tmp(next_state,1:time) =[path_mem(s,1:time-1) i];
                     end;
                 end;
                 s = s+1;
             end;
         end;
         path_metric=path_met_tmp;   % path metric
         path_mem=path_mem_tmp;      % surviving path to each state
     end;
     Data = path_mem(1,1:len_x).';
     clear path_mem;
else            % truncation case
     st=0;
     for time = 2:len_total
         path_met_tmp = inf*ones(n_shift_state,1);
         % calculate branch metrics against every output pattern
         for i = 1:size(out_pat,1);
             branch_metric(i) = sum((out_pat(i,:) - y(time,:)).^2);
         end;
         s=1;
         while s <= n_shift_state
             if path_metric(s) == NaN
                 s = s+1;
             else
                 for i = 0:1     % Input is 0 or 1
                     next_state = transition(s,i+1);
                     branch(next_state) = branch_metric(out_pat_n(s,i+1));
                     if path_met_tmp(next_state) >= path_metric(s)+branch(next_state)
                         path_met_tmp(next_state) = path_metric(s)+branch(next_state);
                         path_mem_tmp(next_state,1:time-st) =[path_mem(s,1:time-st-1) i];
                     end;
                 end;
                 s = s+1;
             end;
         end;
         if time < mem
             path_metric=path_met_tmp;   % path metric
             path_mem=path_mem_tmp;      % surviving path to each state
         else
             path_metric=path_met_tmp;   % path metric
             path_mem=path_mem_tmp;      % surviving path to each state
             smallpath=min(path_mem(:,mem));
             r=1;
             while r <= n_shift_state
                 if path_mem(r,mem) == smallpath
                     x_hat(time-mem+1) = path_mem(r,1);
                     r = n_shift_state;
                 end;
                 r = r+1;
             end;
             path_mem=shift_de(path_mem_tmp,0,1);
             st=st+1;
         end;
     end;
     x_hat = [x_hat path_mem(1,1:mem-1)];
     Data = x_hat(1,1:len_x).';
     clear x_hat;%save memory
     clear path_mem;%save memory
end;