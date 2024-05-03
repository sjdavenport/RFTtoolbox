function sample_vec = vec_ecdf( sample_vec, dist_vec )
% Evaluation of an empirical cumulative distribution function
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  sample_vec   a vector of sample data to find the quantiles of
%  dist_vec     a vector given the data from the distribution that you
%               would like to compare the quantiles of
%--------------------------------------------------------------------------
% OUTPUT
% sample_vec    the ecdf 
%--------------------------------------------------------------------------
% EXAMPLES
% sample_vec = normrnd(0,1,1,10000); dist_vec = normrnd(0,1,1,100000);
% a = vec_ecdf( sample_vec, dist_vec );
% histogram(a)
%
% sample_vec = rlap( 3, [1, 100000]);; dist_vec = rlap( 3, [1, 100000]);;
% a = vec_ecdf( sample_vec, dist_vec );
% histogram(a)
%
% sample_vec = [-1.96, 0, 1.64]; dist_vec = normrnd(0,1,1,100000);
% a = vec_ecdf( sample_vec, dist_vec );
% histogram(a)
%
% sample_vec = [ 1,2]; vec_ecdf( sample_vec, sample_vec+0.1 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------


%%  Main Function Loop
%--------------------------------------------------------------------------

% Sort the vector with the value from the reference distribution
dist_vec = sort(dist_vec);

% Sort the samples from the data and record the sort index
[sample_vec, sort_index] = sort(sample_vec);

% Obtain the inverse sort index so you can put things back in the order
% they came in
sort_index = invPerm(sort_index);

% Calculate the number of samples
nsamples = length(sample_vec);

% Main loop: compare the sample and reference distributions
counter = 1;
for I = 1:nsamples
    found = 0;
    while found == 0
        if sample_vec(I) > dist_vec(counter)
            counter = counter + 1;
            found = 0;
        else
            sample_vec(I) = counter - 1/2;
            found = 1;
        end
        if counter > length(dist_vec)
            break
        end
    end
    % counter
    if counter > length(dist_vec)
        %         sample_vec(I:end) = counter - 1.25;
            sample_vec(I:end) = counter - 1.5;
            % sample_vec(I:end) = counter + 0.5;
            % sample_vec(I:end) = counter - 1;
        break
    end
end

sample_vec = sample_vec/length(dist_vec);
sample_vec = sample_vec(sort_index);

end

