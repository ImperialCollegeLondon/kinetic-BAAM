function [outputs] = evaluate_samples(inputs,parameters,numOut,maxTime,flag1D,flagLbyr)
% A sub-function for evaluating the objective function for all members of a
% population in parallel (with a time limitation per evaluation) using parfeval

N_points = size(inputs, 1); % Number of input points to run [-]
max_duration = maxTime; % Maximum allowed simulation time before termination [s]


outputs = zeros(size(inputs, 1), numOut); % Array for storing model outputs
for j = 1:N_points
    if flagLbyr
    parameters.Lbyr = inputs(j, end);
    inputs(j, 1) = inputs(j, 1) .*parameters.Lbyr./7;
    else
    parameters.rp =    inputs(j, end);
    end
    if flag1D
        F(j) = parfeval(@pvsa_comp, 1, parameters,inputs(j, :));
    else
        F(j) = parfeval(@kBAAM_Outputs_nonIsothermal_dP, 1, parameters,inputs(j, :));
    end
end

completed_sims = zeros(N_points, 1);
all_completed = 0;
tic
while (all_completed == 0) && (toc < max_duration)

    pause(1)

    for j = 1:N_points

        if (convertCharsToStrings(F(1, j).State) == "finished") && ((datenum(F(1, j).FinishDateTime) - datenum(F(1, j).StartDateTime))*86400 < max_duration) && (completed_sims(j) == 0)
            try
                outputs(j, :) = F(1, j).OutputArguments{1, 1}
                completed_sims(j) = 1;
            catch
                cancel(F(1, j))
                completed_sims(j) = 1;
                disp(inputs(j, :))
            end
        end

        if (convertCharsToStrings(F(1, j).State) == "running") && ((now - datenum(F(1, j).StartDateTime))*86400 > max_duration) && (completed_sims(j) == 0)
            cancel(F(1, j))
            completed_sims(j) = 1;
        end

    end

    if completed_sims == ones(N_points, 1)
        all_completed = 1;
    end

end

for j = 1:N_points

    if completed_sims(j) == 0
        cancel(F(1, j))
        fprintf([' fail ', num2str(inputs(j, :))])
    end

end

end