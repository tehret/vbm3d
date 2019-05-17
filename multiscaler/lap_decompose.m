function out = lap_decompose(input, scales, interp, prefix, suffix) 
%
% cellarray = lap_decompose(image, 5)
% cellarray = lap_decompose(image, 5, 'bicubic')
% cellarray = lap_decompose(image, 5, 'bicubic', './pyr', '.tif')
%
if (nargin < 3) || isempty(interp)
     interp = 'lanczos3';
end
if (nargin < 5) 
   savepyr= false;
else
   savepyr = true;
end

out = cell(1,scales);
   for i = 1:scales
      out{i} = input;
      if savepyr==true 
         writeTIFF(out{i}, [prefix, num2str(i-1), suffix]);
      end
      input = imresize(input,0.5,interp);
   end
end


function writeTIFF(data, filename)
% writeTIFF(data, filename)
% writes data as a multi-channel TIFF with single prec. float pixels
   t = Tiff(filename, 'w');
   tagstruct.ImageLength = size(data, 1);
   tagstruct.ImageWidth = size(data, 2);
   tagstruct.Compression = Tiff.Compression.None;
   %tagstruct.Compression = Tiff.Compression.LZW;        % compressed
   tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
   tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
   tagstruct.BitsPerSample =  32;                        % float data
   tagstruct.SamplesPerPixel = size(data,3);
   tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
   t.setTag(tagstruct);
   t.write(single(data));
   t.close();
end
