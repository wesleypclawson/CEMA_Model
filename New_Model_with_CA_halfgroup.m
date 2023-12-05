%%
N = 149;
T = (ceil(N/2));

sizen = [1 2 4 5 8 10 12 15 17 18];
sizen2 = sizen.^2;
noise = 0.5:0.1:1;
n_dec = 0.4:0.1:1;
sl = length(sizen);
ndl = length(n_dec);
nnl = length(noise);

neighbor1 = [1,1,1;1,0,1;1,1,1];
neighbor2 = [1 1 1 1 1; 1 0 0 0 1; 1 0 0 0 1; 1 0 0 0 1; 1 1 1 1 1];
n2wH = 0.25;
n2wC = 0.25;

%cheat_code = 1./(1+exp(-0.025.*(x-324/2)));

R = 1;

figure();
nc = 1;
for ri = 1:R
    for ni = noise(end-2) %*ones(1, length(noise)) %(4)
        results_sz_ndi = zeros([ndl, sl]);
        results_ca = cell([ndl, sl]);

        ndc = 1;
        for ndi = n_dec(4) %*ones(1, length(n_dec)) %(5)
            nec = 1;
            for sz = sizen %(7)
                
                % where the CAs will be held for all runs
                ca = cell(sz^2, R);

                %  Initialize embryos
                for ei = 1:sz^2
                    ca{ei, ri} = zeros(T, N);
                    num_ones = ceil(rho * N);
                    one_index = randsample(N, num_ones, 'false');
                    ca{ei, run}(1, one_index) = 1;
                    %ca{ei, run}(1, 1:num_ones) = 1; < for custom init
                end

                clear('comp_tm1', 'ncomp_tm1')
                comp = ones(sz, sz);
                compt = zeros(sz^2, T);
                compt((end/2)+1:end, :) = 1;
                compt(:, 1) = reshape(comp.',1,[]);

                for t = 1:T
                    if exist('comp_tm1', 'var')
                        comp_tm1 = comp;
                    else
                        comp_tm1 = ones(sz, sz);
                    end

                    % Noise
                    n_hit = rand(size(comp)) < ni;
                    comp(n_hit) = comp(n_hit) - ndi.*comp(n_hit);
                    comp(int8(end/2)+1:end, :) = 1;

                    if t >= 1
                        %  Feedback
                        t2 = zeros(size(comp));
                        for i = 1:sz
                            for j = 1:sz
                                M = zeros(size(comp));
                                M(i, j) = 1;
                                neighbors = comp(conv2(M, neighbor1,'same')>0);
                                neighbors = neighbors(neighbors ~= 1);
                                nneighbors = comp_tm1(conv2(M, neighbor2,'same')>0);
                                nneighbors = nneighbors(nneighbors ~= 1);
                                if comp(i,j) < 1 && comp(i,j) > 0
                                    t2(i, j) = sum([comp(i, j); (1).*neighbors; n2wH.*nneighbors])./ ...
                                        sum([1 (1).*ones(1, length(neighbors)) n2wH.*(ones(1, length(nneighbors)))]);
                                elseif comp(i, j) >= 1
                                    t2(i,j) = 1;
                                elseif comp(i, j) == 0
                                    t2(i, j) = 0;
                                end
                            end
                        end
                        comp = t2;
                    end

                    % Collective Effects
                    if exist('ncomp_tm1', 'var')
                        ncomp_tm1 = ncomp;
                    else
                        ncomp_tm1 = zeros(sz, sz);
                    end

                    if t >= 1
                        t2 = zeros(size(comp));
                        ncomp = ones(size(comp)) - comp;
                        for i = 1:sz
                            for j = 1:sz
                                M = zeros(size(ncomp));
                                M(i, j) = 1;
                                neighbors = ncomp(conv2(M, neighbor1,'same')>0);
                                neighbors = neighbors(neighbors ~= 0);
                                nneighbors = ncomp_tm1(conv2(M, neighbor2,'same')>0);
                                nneighbors = nneighbors(nneighbors ~= 0);
                                if comp(i,j) < 1 && comp(i,j) > 0
                                    t2(i, j) = sum([ncomp(i,j)*ncomp(i,j); ...    % Changed comp to ncomp
                                        ncomp(i,j).*(1).*neighbors; ...
                                        ncomp(i,j).*n2wC.*nneighbors])./ ...
                                        sum([ncomp(i,j) ... %changed comp to ncomp
                                        ncomp(i,j).*(1).*ones(1, length(neighbors)) ...
                                        ncomp(i,j).*n2wC.*(ones(1, length(nneighbors)))]);
                                elseif comp(i, j) >= 1
                                    t2(i,j) = 0;
                                elseif comp(i, j) == 0
                                    t2(i,j) = 1;
                                end
                            end
                        end
                        ncomp = t2;

                        % The component value goes up proportional to it's neighbors
                        % mean distance to one, weighted by its own heath
                        % health value of the component -> comp(comp ~= 1)
                        % distance of each component to 1 -> ncomp(comp ~= 1)
                        % mean distance of neighbors to 1 -> t2(comp ~= 1)
                        %comp(comp ~= 1) = comp(comp ~= 1) + comp(comp ~= 1).*ncomp(comp ~= 1);
                        comp(comp ~= 1) = comp(comp ~= 1) + comp(comp ~= 1).*ncomp(comp ~= 1);
                        comp(comp > 1) = 1;
                    end
                    %comp(comp ~= 1) = comp(comp ~= 1) + ncomp(comp ~= 1);
                    compt(:, t+1) = reshape(comp.',1,[]);

                    % Force the top half of the embryos to always be 1. 
                    if sz > 1
                        comp(int8(end/2)+1:end, :) = 1;
                        %compt((end/2)+1:end, t+1) = 1;
                    end

                    % Update Embryos
                    % Update all embryos
                    for ei = 1:sz^2

                        % Update current CA
                        current_ca = ca{ei, ri};

                        % GKL Rule
                        sr1 = circshift(current_ca(t, :), 1, 2);
                        sr3 = circshift(current_ca(t, :), 3, 2);
                        sl1 = circshift(current_ca(t, :), -1, 2);
                        sl3 = circshift(current_ca(t, :), -3, 2);

                        % If if noise should be applied, do so, if not, update with
                        % GKL
                        if compt(ei, t) <= 0.5
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
                        elseif compt(ei, t) > 0.5
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
                    end
                end
                results_sz_ndi(ndc, nec) = mean(compt(:, T));

                % Go through and save distribution of success values
                success = nan([1, sz^2]);
                for ei = 1:sz^2
                    current_ca = ca{ei, run};
                    success(1, ei) = sum(current_ca(end, :))/(size(current_ca, 2));
                end
                results_ca{ndc, nec} = success;

                nec = nec + 1;
            end
            ndc = ndc + 1;
        end

        nexttile();
        imagesc(results_sz_ndi);
        colorbar;
        caxis([0 1]);
        xticks(1:length(sizen));
        temp = {};
        for i = 1:length(sizen)
            temp{1, end+1} = num2str(sizen(1, i)^2);
        end
        xticklabels(temp);
        yticks(1:length(n_dec));

        temp = {};
        for i = 1:length(n_dec)
            temp{1, end+1} = num2str(n_dec(1, i));
        end
        yticklabels(temp);
        title(['noise = ' num2str(ni)]);
        nc = nc + 1;

    end
end

final = nan(7, 10);
for j = 1:7
    for i = 1:10
        final(j, i) = sum(results_ca{j, i} >= 1)/length(results_ca{j, i});
    end
end

xax = sizen2;
fig = figure(); errorbar(xax, nanmean(final), nanstd(final));
ylim([0 1]);
xlabel('Number of embryos in group');
ylabel('% Non-defect');
fontsize(fig, 10, 'points');
set(gca,'fontname','veranda')

%% 

final = nan(7, 10);
for j = 1:7
    for i = 1:10
        final(j, i) = sum(results_ca{j, i} >= 1)/length(results_ca{j, i});
    end
end

xax = sizen2;
fig = figure(); errorbar(xax, nanmean(final), nanstd(final));
ylim([0 1]);
xlabel('Number of embryos in group');
ylabel('% Non-defect');
fontsize(fig, 10, 'points');
set(gca,'fontname','veranda')