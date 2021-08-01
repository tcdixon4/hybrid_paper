function [unit_data] = add_electrode_depth(unit_data)

%
% Adds a field with the relative depth of each unit, dictated by the
% channel number
%
%
% INPUTS: 
%
% unit_data - struct: 1 x num_units
%             unit-separated data struct containing hemisphere, brain area,
%             firing rate and other metrics
%
% OUTPUTS:
%
% unit_data - struct: 1 x num_units
%             same as unit_data, but with a new field indicating the depth 
%             of the channel that the unit was recorded on and another new 
%             field providing an identifier for the specific probe 
%             insertion that that unit was recorded on (e.g. the left PMd
%             probe on day 1 has it's own identifier)
%
%


%% add depth iteratively

channel2depth = [1:32,...
                 17:32,1:16,...
                 1:32,...
                 17:32,1:16];
channel2probe = [ones(1,32), 2*ones(1,32), 3*ones(1,32), 4*ones(1,32)];
for unit = 1:length(unit_data)
    channel = str2num(unit_data(unit).id(4:6));
    unit_data(unit).depth = channel2depth(channel);
    unit_data(unit).probe_id = ...
        unit_data(unit).session_num*channel2probe(channel);
end


end




