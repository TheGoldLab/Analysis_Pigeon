function dataDir_ = getPigeon_dataDir(name)
% function dataDir_ = getPigeon_dataDir(name)

if nargin < 1 || isempty(name)
    % [~,name] = system('hostname');
    [~,name] = system('ifconfig en0 | grep ether');
end

if strcmpi(name(11:24), '1b:c7:f3:4a:c4') || strcmpi(name(8:24), '4a:86:dd:03:3c:22') || strcmpi(name(8:24), '56:7c:d7:14:e7:56') % strncmpi(name,'Gold',2) || strncmpi(name,'nsc-gold-001.med.upenn.edu',6)
    dataDir_ = fullfile('/Users', 'jigold', 'GoldWorks', 'Mirror_jigold', 'Manuscripts', '2023_Pigeon', 'Data');
else
    fullfile('/Users', 'ishankalburge', 'Documents', 'penn', 'Pigeon', 'Data');
end    