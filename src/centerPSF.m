function [H] = centerPSF(H, thresh)
%CENTERPSF Thresholds and centers nonzero PSF elements within the PSF window
%
% H				a x b x P matrix of PSFs for P images, P=1 for single channel case
% thresh		threshold in the range [0,1] for determining the PSF support (and for actual thresholding, if any). H is constrast-stretched to [0,1] before applying the threshold.
% method		(temporarily removed) thresholding method
%				'none'		(default) PSF is preserved as is, thresholding is used only to determine PSF support
%				'hard'		values < threshold are set to zero
%
% H				a x b x P matrix, each PSF is centered and sums to 1.

hsize = [size(H,1) size(H,2)];

for i=1:size(H,3)
	h = mat2gray(H(:,:,i)); % contrast-stretching to [0,1]	
	
	% nonzero mask
	m = h >= thresh;
	m2 = bwmorph(m, 'clean'); % remove isolated pixels
	if(any(m2(:))) m = m2; end % preserve delta-like PSFs
	
	% determine mask support
	sum1 = sum(m, 1);
	sum2 = sum(m, 2);
	L = [find(sum2, 1, 'first') find(sum1, 1, 'first')];
	R = [find(sum2, 1, 'last') find(sum1, 1, 'last')];
	topleft = fix((L+R+1-hsize)/2); % topleft corner index of the new mask
	
	% indexing (=shifting)
	h = h(max(topleft(1),1):min(topleft(1)+hsize(1)-1,end), max(topleft(2),1):min(topleft(2)+hsize(2)-1,end)); % get the 'existing' data, then pad borders with zeros
	h = padarray(h, max(topleft-[1,1],0), 0, 'post'); % pad with zeros to end up with the same size
	h = padarray(h, max([1,1]-topleft,0), 0, 'pre');
	
	% normalize sum to 1
	H(:,:,i) = h/sum(h(:));
%     H(:,:,i)=h;
end
end

