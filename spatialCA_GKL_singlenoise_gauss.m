rng shuffle
clear

% Set number of cells and time to run
N = 149;
T = (ceil(N/2)+10);

% Set single-run noise
% noise = 0.5;
thresh = 0; % Normal is 0 

for noise = 0.35 %0.4:0.01:0.53 %0.1:0.1:1
    % How many runs (how many repeats do you want to do per noise lvl)
    R = 4;
    e_size = [1 2 3 5 8 10 13 17];
    Embryos = e_size.^2; %[1 5 25 75 100 300];

    % Holders for plotting
    plotmean = nan(1, length(Embryos));
    plotstd = nan(2, length(Embryos));

    for eidx = 1:length(e_size)
        % How many embryos (how many total embroys in the dish)
        num_embryo = e_size(eidx)^2;
        E = num_embryo;

        disp(['Num embryo = ' num2str(num_embryo)]);

        % Spatial distribution
        embryo_id = nan(e_size(eidx), e_size(eidx));
        c = 1;
        for i = 1:e_size(eidx)
            for j = 1:e_size(eidx)
                embryo_id(j, i) = c;
                c = c + 1;
            end
        end

        % Noise and noise tracker
        track_noise = zeros(R, E, T);

        % Set Rho
        rho = 0.6;

        % success ratio
        th_success = 1; %0.9;

        % Main output
        CA_manip = [];
        Survival = zeros(length(noise), R);
        Variability = zeros(length(noise), 1);

        % where the CAs will be held for all runs
        ca = cell(E, R);

        % Track successes over all runs
        success = nan(E, R);

        % Start a progress bar
        pb = 0;
        progressbar('Time to complete');

        % Do R, runs (repetitions of same noise level)
        for run = 1:R

            %  Initialize embryos
            for ei = 1:E
                ca{ei, run} = zeros(T, N);
                num_ones = ceil(rho * N);
                one_index = randsample(N, num_ones, 'false');
                ca{ei, run}(1, one_index) = 1;
                %ca{ei, run}(1, 1:num_ones) = 1; < for custom init
            end

            % track health
            health = ones(e_size(eidx), e_size(eidx), T);

            % Simulate the CAembryo over the defined time period
            for t = 1:T

                % Step 1: See if noise happens
                if t == 1
                    thealth = health(:, :, t);
                    if E == 1
                        rand_poison = rand(e_size(eidx)) <= noise;
                        if rand_poison; thealth = -1; end
                    else
                        rand_poison = rand(e_size(eidx)) <= noise;
                        thealth(rand_poison) = -1;
                    end
                else
                    thealth = health(:, :, t);
                    if E == 1
                        rand_poison = rand(e_size(eidx)) <= noise;
                        if rand_poison; thealth = -1; end
                    else
                        rand_poison = rand(e_size(eidx)) <= noise;
                        thealth(rand_poison) =  thealth(rand_poison)-1*sum(reshape(thealth, [],1') <= thresh)/E;
                        %thealth(rand_poison) = -1*sum(reshape(thealth, [],1') <= thresh)/E;
                    end
                end

                % Step 2: Have the health/disease spread
                %             if t >= 2 && t ~= T
                %                t2health = health(:, :, t-1);
                %                t2health = imgaussfilt(t2health, 2);
                %                health(:, :, t) = health(:, :, t) + t2health; %/eidx;
                %             end
                thealth(thealth > thresh) = thealth(thealth > thresh) + sum(sum(thealth > thresh))./E;
                %thealth(thealth > thresh) = thealth(thealth > thresh) + (e_size(eidx)/13).*sum(sum(thealth > thresh))./E;

                health(:, :, t+1) = imgaussfilt(thealth, 2);

                % Step 3: Update all embryos
                for ei = 1:E

                    % Update current CA
                    current_ca = ca{ei, run};

                    % GKL Rule
                    sr1 = circshift(current_ca(t, :), 1, 2);
                    sr3 = circshift(current_ca(t, :), 3, 2);
                    sl1 = circshift(current_ca(t, :), -1, 2);
                    sl3 = circshift(current_ca(t, :), -3, 2);

                    % If if noise should be applied, do so, if not, update with
                    % GKL
                    %if track_noise(run, ei, t) == 1
                    [a,b] = ind2sub(e_size(eidx), ei);
                    if health(a, b, t) < thresh
                        flipped_bit0 = randsample(3, 1);
                        flipped_bit1 = flipped_bit0;

                        for n = 1:N
                            if current_ca(t, n) == 0
                                update = [sr1(1, n) sr3(1, n) current_ca(t, n)];
                                current_ca(t+1, n) = update(flipped_bit0);
                            elseif current_ca(t, n) == 1
                                update = [sl1(1, n) sl3(1, n) current_ca(t, n)];
                                current_ca(t+1, n) = update(flipped_bit1);
                            end
                        end
                    else
                        for n = 1:N
                            if current_ca(t, n) == 0
                                current_ca(t+1, n) = mode([sr1(1, n) sr3(1, n) current_ca(t, n)]);
                            elseif current_ca(t, n) == 1
                                current_ca(t+1, n) = mode([sl1(1, n) sl3(1, n) current_ca(t, n)]);
                            end
                        end
                    end
                    % End of GKL

                    % Save the updated CA
                    ca{ei, run} = current_ca;

                    if t == T
                        % Check success conditions -> th_success maj
                        if rho > 0.5
                            if sum(current_ca(end, :))/N >= th_success
                                success(ei, run) = 1;
                            else
                                success(ei, run) = 0;
                            end
                        elseif rho < 0.5
                            if sum(current_ca(end, :))/N <= (1-th_success)
                                success(ei, run) = 1;
                            else
                                success(ei, run) = 0;
                            end
                        end
                    end
                    pb = pb + 1;
                    progressbar(pb/(length(noise)*E*R*T));
                end
            end
        end

        Survival(1, :) = sum(success) / E;
        Variability(1, 1) = sum(reshape(success.',1,[])) / (E*R);
        CA_manip.(['n_' num2str(1)]) = ca;

        if num_embryo == 1
            plotmean(1, eidx) = mean(Survival/R,2);
            SEM = std(Survival, [], 2)/length(Survival);
            ts = tinv([0.025 0.975], length(Survival)-1);
            CI = mean(Survival/R,2) + ts*SEM;
            plotstd(:, eidx) = CI';
        else
            plotmean(1, eidx) = mean(Survival,2);
            SEM = std(Survival, [], 2)/length(Survival);
            ts = tinv([0.025 0.975], length(Survival)-1);
            CI = mean(Survival,2) + ts*SEM;
            plotstd(:, eidx) = CI';
        end
    end

    errorbar(Embryos, plotmean, plotstd(1, :)/10, plotstd(2, :)/10); hold on;
end
%%

legendCell = cellstr(num2str(num_embryo', 'N=%-d'));
legend(legendCell{:});


