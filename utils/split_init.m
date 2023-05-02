function init = split_init(alfarange, guess, nx, nu)
    init.index = [];
    for i = 1:length(alfarange)
        [~, index] = min(abs(guess.q(1:nx:end) - alfarange(i)));
        init.index = [init.index;index];
        init.q(i,:) = guess.q(nx*(index-1)+1:nx*(index));
        if i > 1
            if index > 1
                rw = guess.w((init.index(i-1)-1)+1:(index-1),:);
                init.w(i-1,:) = mean(rw,1);
            else
                init.w(i-1,:) = guess.w(index,:);    
            end

        end
    end
end