function dataDir_ = getPigeon_dataDir(name)
% function dataDir_ = getPigeon_dataDir(name)

if nargin < 1 || isempty(name)
    [~,name] = system('hostname');
end

if strncmpi(name,'Gold',2) || strncmpi(name,'nsc-gold-001.med.upenn.edu',6)
    dataDir_ = fullfile('/Users', 'jigold', 'GoldWorks', 'Mirror_jigold', 'Manuscripts', '2023_Pigeon', 'Data');
else
    fullfile('/Users', 'ishankalburge', 'Documents', 'penn', 'Pigeon', 'Data');
end    