using Images, TestImages, ImageView, FileIO

function HorizontalRF(I, D, sigma)
	
	a = exp(-sqrt(2) / sigma)
	
	F = deepcopy(I)
	
	V = a .^ D
	
	(numberOfChannels, height, width) = size(I)
	
	for i = 2:width
		for c = 1:numberOfChannels
			F[c,:,i] = F[c,:,i] .+ V[:,i] .* ( F[c,:,i-1] .- F[c,:,i] )
		end
	end
	
	for i = width-1:-1:1
		for c = 1:numberOfChannels
        	F[c,:,i] = F[c,:,i] .+ V[:,i+1] .* ( F[c,:,i+1] .- F[c,:,i] );
		end
	end
	
	return F
	
end

function image_transpose(I)
	(num_channels, h, w) = size(I);
    
    T = zeros(num_channels, w, h);
    
    for c = 1:num_channels
        T[c,:,:] = I[c,:,:]';
	end
	
	return T
	
end


img = testimage("mandrill")

img = RGB{Float64}.(img)

sigma_s = 60;
sigma_r = 0.4;
num_iterations = 3

(h,w) = size(img)

dyPartial = diff(img,dims=1)
dxPartial = diff(img,dims=2)
dxPartial = channelview(dxPartial)
dyPartial = channelview(dyPartial)
dIdx = zeros(h,w)
dIdy = zeros(h,w)

for c = 1:3
	for y = 1 : h-1
	  for x = 1 : w-1
		   dIdx[y,x+1] = dIdx[y,x+1] + abs(dxPartial[c,y,x]);
		   dIdy[y+1,x] = dIdy[y+1,x] + abs(dyPartial[c,y,x]);
		end
	end
end

#Compute the derivatives of the horizontal and vertical domain transforms.
dHdx = (1 .+ sigma_s/sigma_r * dIdx);
dVdy = (1 .+ sigma_s/sigma_r * dIdy);
	
dVdy = dVdy';

sigma_H = sigma_s;

F = deepcopy(channelview(img))
N = num_iterations

for i = 0:num_iterations-1

	  #Compute the sigma value for this iteration
	  sigma_H_i = sigma_H * sqrt(3) * 2^(N - (i + 1)) / sqrt(4^N - 1);

	  global F = HorizontalRF(F, dHdx, sigma_H_i);
	  global F = image_transpose(F);
	  global F = HorizontalRF(F, dVdy, sigma_H_i);
	  global F = image_transpose(F);

end

F = colorview(RGB, F)

imshow(F)


