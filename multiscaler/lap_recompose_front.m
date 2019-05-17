function ret = lap_recompose_front(prefix, levels, suffix, g, tau, cur, interp)
%
%  x = lap_recompose_front(filename, levels, '.tif', recfactor,0,0);
%  x = lap_recompose_front(filename, levels, '.tif', recfactor,0,0, 'lanczos3');
%
if (nargin < 7) || isempty(interp)
     interp = 'lanczos3';
end

pyr={};
for cur=0:levels-1
   name = [prefix,num2str(cur),suffix];
   t = Tiff(name);
   pyr{cur+1} = double(t.read());
end

ret = lap_recompose(pyr, g, tau, interp);

end
