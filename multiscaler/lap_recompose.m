function ret = lap_recompose(pyr, g, tau, interp)
%
%  x = lap_recompose(pyr, 0.8, 0.0);
%  x = lap_recompose(pyr, 0.8, 0.0, 'lanczos3');
%
if (nargin < 4) || isempty(interp)
     interp = 'lanczos3';
end

image = pyr{1};

if (numel(pyr)==1) % if last level
   ret=image;
   return 
end

yL = lap_recompose({pyr{2:end}}, g, tau, interp) ;

yH = image;

% pyramid recomposition
sz=size(yH);
H = yH - up(gblur(down(yH,interp),g),sz,interp);
L = up(gblur(yL,g),sz,interp);

% threshold high frequency values larger than tau
H = (abs(H) >= tau) .* H;
ret  = L + H;
%writeTIFF(ret, [prefix,num2str(cur),'x',suffix]);

end


function ret=up(image,sz,interp)
   if (nargin < 3) || isempty(interp)
        interp = 'lanczos3';
   end
   ret = imresize(image,2,interp);
   sz2=size(ret);

   if sz2(1)>sz(1) 
       ret = ret(1:end-1,:,:);
   end
   if sz2(1)<sz(1) 
       ret = ret([1:end,end],:,:);
   end
   if sz2(2)<sz(2) 
       ret = ret(:,[1:end,end],:);
   end
   if sz2(2)>sz(2) 
       ret = ret(:,1:end-1,:);
   end
end

function ret = down(image,interp)
   if (nargin < 2) || isempty(interp)
        interp = 'lanczos3';
   end
   ret = imresize(image,0.5,interp);
end

function R = gblur(I, s)
    if s == 0 
        R = I;
    else
    border_mode = 'symmetric';
    supp = max(floor(s)*2,5);
    filter = fspecial( 'gaussian', supp,s);
    % low pass, convolve with separable filter
    R = imfilter(I,filter,border_mode);     %horizontal
    R = imfilter(R,filter',border_mode);    %vertical
  % R = imgaussfilt(I,s,'Padding','symmetric');
    end
end
