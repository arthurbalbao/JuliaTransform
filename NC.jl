using Images, TestImages, ImageView, FileIO

img = load("statue.png")

img = RGB{Float32}.(img)


#extracts partial derivates from image
diff(img,1)
