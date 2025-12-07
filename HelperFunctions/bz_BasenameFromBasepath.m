function basename = bz_BasenameFromBasepath(basepath)
    % Extract the basename from the given basepath
    [~, basename, ~] = fileparts(basepath);
end
