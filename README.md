# About
This project uses local warping to make panoramas with irregular boundaries rectangular. 
It is based on the algorithm described in ["Rectangling Panoramic Images via Warping"](http://research.microsoft.com/en-us/um/people/kahe/publications/sig13pano.pdf) by Kaiming He, Huiwen Chang, and Jian Sun.

## Algorithm
The algorithm uses seam insertion to locally warp the image to a rectangle, imposes a displacement mesh over the image, calculates an equation representing the "energy" of the mesh, optimizes the energy function to yield a final mesh, and produces a rectangular output image by interpolating from this mesh.

## Implementation
This implementation works, but its energy function is missing one term described in the paper, line preservation. 
This energy function preserves shapes and fits the mesh to a rectangle, but it does not preserve lines. 
A significant portion of the code to calculate line preservation is included but is not currently being used.

# How to use
This application requires the input image to be a simple 24-channel bitmap, and the unfilled space around the image must be white (#FFFFFF).
Run the application from XCode, right-click in the application to open an image, and select "Rectangle". You can then save the output as a bitmap.
