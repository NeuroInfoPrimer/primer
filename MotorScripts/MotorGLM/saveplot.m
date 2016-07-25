function saveplot(h, filename, format, dimensions)
%saveplot	Saves the specified plot.
%
% Usage:
%			saveplot(h, filename, format, dimensions, fontsize)
%
% Input:
%			h = figure object. Can just use 'gcf'
%			filename = output filename
%			format = one of 'eps', 'png', 'jpg', 'pdf' or 'svg' (optional, default 'eps')
%			dimensions = (optional, default = [6, 4]) width and height of figure in inches
%
% Examples:
%			plot(0:10, (0:10).^2);
%			xlabel('time');
%			saveplot(gcf, './test.eps', 'eps');

%Note that matlab saveas and open can be used to store figures for editing later...

	if (nargin < 3)
		format = 'eps';
		dimensions = [6 4];
	elseif nargin < 4
		dimensions = [6 4];
	elseif nargin == 4
		if length(dimensions) == 1
			throw(MException('Argument:Error', 'Must specify a tuple of image dimension'));
		end
	end

	pts = 4*1/72;
	box = [pts, pts, dimensions(1)-2*pts, dimensions(2)-2*pts];

	renderer = '-painters';
	if strcmp(format, 'eps')
		dev = '-depsc';
	elseif strcmp(format, 'png')
		dev = '-dpng';
		renderer = '-zbuffer';
	elseif strcmp(format, 'svg')
		dev = '-dsvg';
		%renderer = '-zbuffer';
	elseif strcmp(format, 'pdf')
		dev = '-dpdf';
	elseif strcmp(format, 'jpg')
		dev = '-djpeg';
		%renderer = '-opengl';
	else
		throw(MException('ArgumentError','Format must be one of eps or png or pdf or jpeg'));
	end

	set(h, 'paperunits', 'inches');
	set(h, 'papersize', dimensions);
	set(h, 'paperposition', box);
	set(h, 'paperpositionmode', 'manual');
	print(h, dev, '-r300', renderer, filename);
end
