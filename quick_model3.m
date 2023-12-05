%% Plot stuff
colors = [
    255,247,188
    254,227,145
    254,196,79
    254,153,41
    236,112,20
    204,76,2
    153,52,4
    102,37,6
    0, 0, 0]./255;
colors = flip(colors);

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

figure();
nc = 1;
for ni = noise(2)
    results_sz_ndi = zeros([ndl, sl]);

    ndc = 1;
    for ndi = n_dec(5)
        nec = 1;
        for sz = sizen(9)
            % establish weights based on size
            
%             n2wH = 1-(1./(1+exp(-0.025.*((sz^2)-324/2))));
%             n2wC = 1./(1+exp(-0.025.*((sz^2)-324/2)));
            
            clear('comp_tm1', 'ncomp_tm1')
            comp = ones(sz, sz);
            compt = zeros(sz^2, T);
            compt(:, 1) = reshape(comp.',1,[]);

            for t = 1:T
                if exist('comp_tm1', 'var')
                    comp_tm1 = comp;
                else
                    comp_tm1 = ones(sz, sz);
                end

%                 n2wH = sum(compt(:, t) < 1 & compt(:, t) > 0) / sz^2; 
%                 n2wC = n2wH/2;

                % Noise
                n_hit = rand(size(comp)) < ni;
                comp(n_hit) = comp(n_hit) - ndi.*comp(n_hit); 

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
%                     for i = 1:sz
%                         for j = 1:sz
%                             M = zeros(size(ncomp));
%                             M(i, j) = 1;
%                             neighbors = ncomp(conv2(M, neighbor1,'same')>0);
%                             neighbors = neighbors(neighbors ~= 0);
%                             nneighbors = ncomp_tm1(conv2(M, neighbor2,'same')>0);
%                             nneighbors = nneighbors(nneighbors ~= 0);
%                             if comp(i,j) < 1 && comp(i,j) > 0
%                                 t2(i, j) = sum([ncomp(i,j)*ncomp(i,j); ...    % Changed comp to ncomp
%                                                ncomp(i,j).*(1).*neighbors; ...
%                                                ncomp(i,j).*n2wC.*nneighbors])./ ...
%                                     sum([ncomp(i,j) ... %changed comp to ncomp
%                                     ncomp(i,j).*(1).*ones(1, length(neighbors)) ...
%                                     ncomp(i,j).*n2wC.*(ones(1, length(nneighbors)))]);
%                             elseif comp(i, j) >= 1
%                                 t2(i,j) = 0;
%                             elseif comp(i, j) == 0
%                                 t2(i,j) = 1;
%                             end
%                         end
%                     end
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

            end

            results_sz_ndi(ndc, nec) = mean(compt(:, T));
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
