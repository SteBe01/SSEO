function [data_med, data_sig, u] = data4(dataM)
%DATA4 computes the medium value, standard deviation and useful percentage
%of the values collected in dataM (4-Dimensional tensor input)
%
%INPUT:
% dataM : [r x c x p x s] 4-D tensor collecting the datas (any dimension)
%
%OUTPUTS:
% data_med : [-] medium value of the useful datas (does not include NaN in the computation
% data_sig : [-] standard deviation of the datas 
% u        : [-] fraction of useful datas (number of non NaN numbers over (r*c*p*s) )
%
% outputs dimensions are consistent with the input dimension
%
% if one or more datas are equal to Inf, the medium value is Inf, as its
% standard deviation.
%
%AUTHORS: Ancillotti G., Tartaglia D., Tessarollo A., Bolsi P.

[r, c, p, s]=size(dataM);

data_sum=0;
x=0;
for i=1:r
    for j=1:c
        for k=1:p
            for t=1:s
                if 1==isnan( dataM(i, j, k, t) )

                else
                    data_sum=data_sum+dataM(i, j, k, t);
                    x=x+1;
                end
            end
        end
    end
end

data_med=data_sum/x;
u=x/(r*c*p*s);

data_sig2=zeros(x, 1);
y=0;

for i=1:r
    for j=1:c
        for k=1:p
            for t=1:s
                if 1==isnan( dataM(i, j, k, t) )

                else
                    y=y+1;
                    data_sig2(y)=(dataM(i, j, k, t) - data_med)^2;
                end
            end
        end
    end
end

data_sig2_med= sum(data_sig2)/x;

data_sig= sqrt( data_sig2_med );

end

