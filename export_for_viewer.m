function export_for_viewer(D, leaders, filename)
%EXPORT_FOR_VIEWER  Save a deployment and leader choices for the Swarm Rally Race web viewer.
%
%   export_for_viewer(D, S, filename)
m = 4;
seed = 7;                        
D = make_deployment(seed);
L(1).name = 'Rohan_Jordan_Siddharth Strategy';
L(1).S    = choose_leaders(D.G, D.P0, D.faulty, D.g, m);

L(2).name = 'Random';
L(2).S    = leaders_random(D.G, D.P0, D.faulty, D.g, m);

L(3).name = 'Highest degree';
L(3).S    = leaders_degree(D.G, D.P0, D.faulty, D.g, m);

L(4).name = 'Closest to g';
L(4).S    = leaders_closest(D.G, D.P0, D.faulty, D.g, m);

export_for_viewer(D, L, ['seed' num2str(seed) '.json'])
%
%   Then open the viewer page and load the file. The page re-runs the
%   Part III model (lambda = 0.5 for leaders, T = 150) and shows each
%   strategy's score J, which should match your own score_mission.

if nargin < 3
    filename = 'swarm.json';
end
if isnumeric(leaders)
    leaders = struct('name', 'My leaders', 'S', leaders);
end
assert(isstruct(leaders) && isfield(leaders, 'name') && isfield(leaders, 'S'), ...
    'leaders must be a vector of indices or a struct array with fields name and S.');
assert(numel(leaders) <= 6, 'The viewer shows at most 6 strategies at once.');

[i, j] = find(triu(D.G));
out.format  = 'ese105-swarm-v1';
out.seed    = D.seed;
out.r       = D.r;
out.P0      = D.P0;
out.edges   = [i j];
out.faulty  = D.faulty(:)';
out.g       = D.g(:)';
out.leaders = struct('name', {leaders.name}, 'S', cellfun(@(s) s(:)', {leaders.S}, 'UniformOutput', false));

fid = fopen(filename, 'w');
assert(fid > 0, 'Could not open %s for writing.', filename);
fprintf(fid, '%s', jsonencode(out));
fclose(fid);
fprintf('Saved %s. Load it in the Swarm Rally Race viewer.\n', filename);
end
