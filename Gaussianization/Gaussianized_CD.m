nvals = 100000;
mu = 0:0.02:1;
X = normrnd(0,1,length(mu),nvals);

f = @(x) (abs(x).^(1/4)).*sign(x);

mX = mean(f(X+ mu'), 2);
sX = std(f(X+mu'), 0, 2);

plot(mu, mX./sX)
hold on
plot(mu, mu)
legend('ratio', 'mean', 'Location', 'NW')

%%
nvals = 100000;
mu = 0:0.1:1;
X = normrnd(0,1,length(mu),nvals);

f = @(x) (abs(x).^(3/4)).*sign(x);

mX = mean(f(X+ mu'), 2);
sX = std(f(X+mu'), 0, 2);

plot(mu, mX./sX)
hold on
plot(mu, mu)
legend('ratio', 'mean', 'Location', 'NW')

%% Other functions
nvals = 100000;
mu = 0:0.02:1;
X = normrnd(0,1,length(mu),nvals); 

f = @(x) norminv(tcdf(x, 3));
f = @(x) sign(x).*abs(norminv(tcdf(x, 3))).^(3/4);

mX = mean(f(X+ mu'), 2);
sX = std(f(X+mu'), 0, 2);

plot(mu, mX./sX)
hold on
plot(mu, mu)
legend('ratio', 'mean', 'Location', 'NW')

%% t random variables
df = 3;
nvals = 10000;
mu = 0:0.02:1;
% mu = 0.:0.001:0.01;
% X = normrnd(0,1,length(mu),nvals);
X = trnd(df,length(mu),nvals);

% f = @(x) norminv(tcdf(x, df));
f = @(x) (abs(x).^(1/2)).*sign(x);

mX = mean(f(X+ mu'), 2);
sX = std(f(X+mu'), 0, 2);

plot(mu/mean(std(X,0,2)), mX./sX)
hold on
plot(mu, mu)
legend('ratio', 'mean', 'Location', 'NW')


%%
plot(mu, mX./sX - mu')
%%
x = -10:0.01:10;
plot(x,f(x))
hold on 
plot(x,x,'--')