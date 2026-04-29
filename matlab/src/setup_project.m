function paths = setup_project()
%SETUP_PROJECT Add local paths and create assignment output folders.

srcDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(srcDir);
rootDir = fileparts(matlabDir);

paths = struct();
paths.root = rootDir;
paths.matlab = matlabDir;
paths.src = srcDir;
paths.figures = fullfile(matlabDir, 'figures');
paths.results = fullfile(matlabDir, 'results');
paths.structopt = fullfile(rootDir, 'structopt');
paths.lectureChapter9 = fullfile(rootDir, 'structopt_exercises_lectures', 'chapter 9');
paths.stabil = fullfile(rootDir, 'stabil', 'stabil-3.1');

addpath(paths.src);

if isfolder(paths.structopt)
  addpath(paths.structopt);
end

if isfolder(paths.lectureChapter9)
  addpath(paths.lectureChapter9);
end

if ~isfolder(paths.figures)
  mkdir(paths.figures);
end

if ~isfolder(paths.results)
  mkdir(paths.results);
end

end

