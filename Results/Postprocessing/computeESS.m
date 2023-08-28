function Kess = computeESS(acf)

K = length(acf)-1;
jump = 1000;
for i = 1:1000:K
    if mod(i,10000) == 0
        disp(['10000/',num2str(K), ' done']);
    end
    tau = 3*(2*sum(acf(1:(i+1)))-1);
    if i >= tau
        break;
    end
end

for j = (i-jump):i
    tau = 3*(2*sum(acf(1:(j+1)))-1);
    if j >= tau
        break;
    end
end

Kess = K/tau;