function [stimulus] = ReadFramev2(fid,Nsamp,Nd,Nv,cx,x0,y0)

% Parameters:   fid - the file identifier of the .raw file
%               Nsamp - number of frames in the film
%               Nd - pixel length of one side of a square frame
%               Nv - pixel length of one side of patched square stimulus
%               cx - spatial downsampling; 1 is no downsampling, 2 averages
%                   2 pixels, etc.
%               x0 - x coordinate of the corner of the patch in the full
%                   stimulus space closest to the y-axis
%               y0 - y coordinate of the corner of the patch in the full
%                   stimulus space closest to the x-axis

fsize = Nv^2/cx^2;
stimulus = zeros(Nsamp*fsize,1);

for j=0:Nsamp-1
    stimulus_full_frame = fread(fid,Nd^2,'uint8');
    stimulus_full_frame = reshape(stimulus_full_frame,Nd,Nd)';
    patch_final = stimulus_full_frame(y0:y0+Nv-1,x0:x0+Nv-1);
    stimulus(j*fsize+1:(j+1)*fsize) = reshape(patch_final,fsize,1);
end