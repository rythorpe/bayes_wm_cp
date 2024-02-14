function estimate=fixedLRmodel(outcomes, learningRate, initialGuess)
% simple fixed learning rate delta rule updating. vanilla.


estimate=nan(size(outcomes));
estimate(1)=initialGuess;
for i = 1:length(outcomes)-1
    PE=outcomes(i)-estimate(i);
    estimate(i+1)=estimate(i)+PE.*learningRate;
end


