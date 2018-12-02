using Images, TestImages, ImageView, FileIO

# same code as NC filter

function HorizontalRF(Image, Derivative, sigma)
	
	a = exp(-sqrt(2) / sigma)
	
	F = deepcopy(Image)
	
	V = a .^ Derivative
	
	(numberOfChannels, height, width) = size(Image)
	
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