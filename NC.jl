using Images, TestImages, ImageView, FileIO

function F = TransformedDomainBoxFilter_Horizontal(I, xform_domain_position, box_radius)

  I = channelview(I)

  (h w channels) = size(I)
  h


end

function NC(sigma_s, sigma_r, num_iterations)

sigma_s = 60;
sigma_r = 0.4;
num_iterations = 3
img = load("statue.png")

img = RGB{Float64}.(img)

(h,w) = size(img)

#extracts partial derivates from image

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

  #Integrate the domain transforms.
  ct_H = cumsum(dHdx, dims=3);
  ct_V = cumsum(dVdy, dims=2);

  #transpose ct vertical
  ct_V = ct_V'

   sigma_H = sigma_s;

    i = 0
    N = num_iterations

       #Compute the sigma value for this iteration
      sigma_H_i = sigma_H * sqrt(3) * 2^(N - (i + 1)) / sqrt(4^N - 1);

      #Compute the radius of the box filter with the desired variance.
      box_radius = sqrt(3) * sigma_H_i;

#      F = TransformedDomainBoxFilter_Horizontal(F, ct_H, box_radius);
#      F = image_transpose(F);

      # Compute the lower and upper limits of the box kernel in the transformed domain.
         l_pos = ct_H .- box_radius;
         u_pos = ct_H .+ box_radius;

         # Find the indices of the pixels associated with the lower and upper limits
         # of the box kernel.


         l_idx = zeros(size(ct_H));
         u_idx = zeros(size(ct_H));

         for row = 1:h
             xform_domain_pos_row = ct_H[row,:];

             l_pos_row = l_pos[row,:];
             u_pos_row = u_pos[row,:];

             local_l_idx = zeros(1, w);
             local_u_idx = zeros(1, w);

             local_l_idx[1] = findfirst(x -> x>l_pos_row[1], xform_domain_pos_row)
             local_u_idx[1] = findfirst(x -> x>u_pos_row[1], xform_domain_pos_row)

             for col = 2:w
                 nextLowIndex = findfirst(x -> x>l_pos_row[col], xform_domain_pos_row[(Int.(local_l_idx[col-1]):end)])
                 nextHighIndex = findfirst(x -> x>u_pos_row[col], xform_domain_pos_row[(Int.(local_u_idx[col-1]):end)])
                 local_l_idx[col] = local_l_idx[col-1] + nextLowIndex - 1
                 local_u_idx[col] = local_u_idx[col-1] + nextHighIndex - 1
             end

             l_idx[row,:] = local_l_idx;
             u_idx[row,:] = local_u_idx;
         end

         # Compute the box filter using a summed area table.
         SAT            = zeros(3,h,w+1);
         SAT[:,:,2:end] = cumsum(I, dims=3);
         row_indices    = repeat((1:h), 1, w)
         F              = zeros(3,h,w);

         for c = 1:3
             a = sub2ind(size(SAT), row_indices, l_idx, repmat(c, h, w));
             b = sub2ind(size(SAT), row_indices, u_idx, repmat(c, h, w));
             F[c,:,:] = (SAT[b] - SAT[a]) ./ (u_idx - l_idx);
         end



     # F = TransformedDomainBoxFilter_Horizontal(F, ct_V, box_radius);
      #F = image_transpose(F);
      end
