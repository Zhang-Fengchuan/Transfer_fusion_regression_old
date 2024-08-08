function S = mcp_function(t,a,lambda)
if t >= a*lambda
    S = 0.5*a*lambda.^2;
else
    S = lambda*t-0.5*a*lambda.^2;
end
end