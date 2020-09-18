X = normrnd(0,1,10,nvals);

f = @(x) norminv(tcdf(x, 3));

mu = 0;
subplot(2,1,1)
histogram(X+mu);
subplot(2,1,2)
histogram(f(X+mu))

%%
df = 10;
X = trnd(df,10,nvals);

f = @(x) norminv(tcdf(x, df));

mu = 0.5;
subplot(2,1,1)
histogram(X+mu);
subplot(2,1,2)
histogram(f(X+mu))