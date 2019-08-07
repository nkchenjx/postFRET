upper = zeros(config.numState);
lower = zeros(config.numState);

for i = 1:config.numState
    for j = 1:config.numState
        a = rateHistory(i,j,:);
        a = a(:);
        m = mean(a);
        b = a(a>m); %upper
        c = a(a<m); %lower
        upper(i,j) = sqrt(sum((b-m).^2)/length(b));
        lower(i,j) = sqrt(sum((c-m).^2)/length(c));
    end
end

upper
lower