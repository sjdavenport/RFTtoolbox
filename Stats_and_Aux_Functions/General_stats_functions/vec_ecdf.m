function sample_vec = vec_ecdf( sample_vec, dist_vec )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  sample_vec   a vector of sample data to fnd the quantiles of
%  dist_vec     a vector given the data from the distribution that you
%               would like to compare the quantiles of
%--------------------------------------------------------------------------
% OUTPUT
% sample_vec
%--------------------------------------------------------------------------
% EXAMPLES
% sample_vec = normrnd(0,1,1,10000); dist_vec = normrnd(0,1,1,100000);
% a = vec_ecdf( sample_vec, dist_vec );
% histogram(a)
%
% sample_vec = [-1.96, 0, 1.64]; dist_vec = normrnd(0,1,1,100000);
% a = vec_ecdf( sample_vec, dist_vec );
%
% sample_vec = [ 1,2]; vec_ecdf( sample_vec, sample_vec+0.1 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------


%%  Main Function Loop
%--------------------------------------------------------------------------
dist_vec = sort(dist_vec);
[sample_vec, sort_index] = sort(sample_vec);
sort_index = invPerm(sort_index);

nsamples = length(sample_vec);

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
    if counter > length(dist_vec)
        %         sample_vec(I:end) = counter - 1.25;
            sample_vec(I:end) = counter - 1.5;
        break
    end
end

sample_vec = sample_vec/length(dist_vec);
sample_vec = sample_vec(sort_index);

end

