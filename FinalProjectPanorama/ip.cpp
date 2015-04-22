//
//  ip.cpp
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/20/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define INFT         1e8

/*
 * Applies the given effect to each pixel in the given source image.
 */
Image* ip_apply_effect_per_pixel (Image* src, Pixel (*effect)(Pixel, void*), void* args)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            Pixel outputPixel = (*effect)(srcPixel, args);
            output->setPixel(x, y, outputPixel);
        }
    }
    return output;
}


/*
 * Generates a solid-color image with the color of the given pixel.
 */
Image* ip_generate_solid_image (int width, int height, Pixel pixel)
{
    Image* output = new Image(width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            output->setPixel(x, y, pixel);
        }
    }
    return output;
}

/*
 * Calculates the psychosomatic intensity.
 */
double ip_get_intensity (Pixel srcPixel)
{
    const double redIntensityWeight = 0.2126;
    const double greenIntensityWeight = 0.7152;
    const double blueIntensityWeight = 0.0722;
    
    double red = srcPixel.getColor(RED);
    double green = srcPixel.getColor(GREEN);
    double blue = srcPixel.getColor(BLUE);
    return clamp(redIntensityWeight * red +
                 greenIntensityWeight * green +
                 blueIntensityWeight*blue,
                 0, 1);
}

/*
 * Creates a pixel with as gray of its psychosomatic intensity.
 */
Pixel ip_gray_pixel (Pixel srcPixel, void* args)
{
    double intensity = ip_get_intensity(srcPixel);
    Pixel outputPixel = Pixel(intensity, intensity, intensity);
    return outputPixel;
}


/*
 * Creates a new image with values equal to the psychosomatic intensities
 * of the source image.
 */
Image* ip_gray (Image* src)
{
    return ip_apply_effect_per_pixel(src, &ip_gray_pixel, NULL);
}

/*
 * Transposes the image between normal and 90 degrees CW and mirror (so top is top or left is top).
 * The mirroring is to preserve +x => +y.
 */
Image* ip_transpose(Image* src, bool forward)
{
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    Image* output = new Image (srcHeight, srcWidth);
    if (forward) {
        for (int x = 0; x < srcHeight; x++) {
            for (int y = 0; y < srcWidth; y++) {
                Pixel srcPixel = src->getPixel(y, srcHeight-1-(srcHeight - 1 - x));
                output->setPixel(x, y, srcPixel);
            }
        }
    } else {
        for (int x = 0; x < srcHeight; x++) {
            for (int y = 0; y < srcWidth; y++) {
                Pixel srcPixel = src->getPixel(srcWidth-1-(srcWidth - 1 - y), x);
                output->setPixel(x, y, srcPixel);
            }
        }
    }
    return output;
}

/*
 * Returns the energy function of pixel at the given coordinate in the given image,
 * defined as the sum of the change in intensity from the pixel's x-neighbor plus the change
 * in intensity from the pixel's y-neighbor (using x+1 and y+1 whenever possible, else x-1 and y-1).
 * The double returned will be in the range [0, 2].
 */
double ip_get_energy(Image* src, int x, int y)
{
    int width = src->getWidth();
    int height = src->getHeight();
    double energy = 0;
    Pixel srcPixel = src->getPixel(x, y);
    if (x + 1 < width) {
        Pixel rightNeighbor = src->getPixel(x + 1, y);
        energy += fabs(ip_get_intensity(rightNeighbor) - ip_get_intensity(srcPixel));
    } else if (x - 1 >= 0) {
        Pixel leftNeighbor = src->getPixel(x - 1, y);
        energy += fabs(ip_get_intensity(srcPixel) - ip_get_intensity(leftNeighbor));
    }
    if (y + 1 < height) {
        Pixel bottomNeighbor = src->getPixel(x, y + 1);
        energy += fabs(ip_get_intensity(bottomNeighbor) - ip_get_intensity(srcPixel));
    } else if (y - 1 >= 0) {
        Pixel topNeighbor = src->getPixel(x, y - 1);
        energy += fabs(ip_get_intensity(srcPixel) - ip_get_intensity(topNeighbor));
    }
    return energy;
}

/*
 * Returns the forward energy given a pixel and its neighbor that are being removed.
 * This method only deals with vertical seams, since it expects horizontal images to be flipped.
 */
double ip_get_forward_energy(Image* src, int x, int y, NeighborDirection direction)
{
    int width = src->getWidth();
    double energy = 0;
    if (direction == DIRECTION_TOP_LEFT) {
        if (x <= 0 || y <= 0) {
            return energy;
        }
        Pixel bottomLeftPixel = src->getPixel(x - 1, y);
        // add top and right
        Pixel topPixel = src->getPixel(x, y - 1);
        energy += fabs(ip_get_intensity(topPixel) - ip_get_intensity(bottomLeftPixel));
        if (x + 1 < width) {
            Pixel rightPixel = src->getPixel(x + 1, y);
            energy += fabs(ip_get_intensity(rightPixel) - ip_get_intensity(bottomLeftPixel));
        }
        return energy;
    } else if (direction == DIRECTION_TOP) {
        // add left and right
        if (x - 1 >= 0 && x + 1 < width) {
            Pixel leftPixel = src->getPixel(x - 1, y);
            Pixel rightPixel = src->getPixel(x + 1, y);
            energy += fabs(ip_get_intensity(leftPixel) - ip_get_intensity(rightPixel));
        }
        return energy;
    } else if (direction == DIRECTION_TOP_RIGHT) {
        if (x >= width - 1 || y <= 0) {
            return energy;
        }
        Pixel bottomRightPixel = src->getPixel(x + 1, y);
        // add top and left
        Pixel topPixel = src->getPixel(x, y - 1);
        energy += fabs(ip_get_intensity(topPixel) - ip_get_intensity(bottomRightPixel));
        if (x - 1 >= 0) {
            Pixel leftPixel = src->getPixel(x - 1, y);
            energy += fabs(ip_get_intensity(leftPixel) - ip_get_intensity(bottomRightPixel));
        }
        return energy;
    } else {
        // this case shouldn't happen
        return 0;
    }
}

/*
 * Shows the energy of each pixel.
 */
Image* ip_energy(Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            double energy = ip_get_energy(src, x, y);
            Pixel outputPixel = Pixel(energy / 2, energy / 2, energy / 2);
            output->setPixel(x, y, outputPixel);
        }
    }
    return output;
}

/*
 * Determines whether the given pixel is 'transparent'.
 */
bool ip_is_transparent(Pixel pixel) {
    return (pixel.getColor(RED) == 1 && pixel.getColor(GREEN) == 1 && pixel.getColor(BLUE) == 1);
}

/*
 * Determines whether the given pixel is part of the 'transparent' border.
 */
bool ip_is_border(Image* src, int x, int y) {
    Pixel pixel = src->getPixel(x, y);
    if (!ip_is_transparent(pixel)) {
        return false;
    }
    
    
    1;return true; //TODO: BFS
    // do a BFS to reach a border pixel
    
//    expandedPixels = set()
//    fringe = util.Queue()
//    fringe.push(startPos)
//    while not fringe.isEmpty():
//        offPixel = fringe.pop()
//        if offPixel not in expandedPixels:
//            expandedPixels.add(offPixel)
//            for neighbor in getOffNeighbors(datum, offPixel):
//                fringe.push(neighbor)
//                return expandedPixels
}

/*
 * Computes the lowest-energy seam in the range [start, end] using forward energy and returns the
 * coordinates of its pixels. Start is top/left, end is bottom/right.
 */
int* ip_get_seam_in_range(Image* src, SeamOrientation orientation, int start, int end)
{
    if (orientation == ORIENTATION_HORIZONTAL) {
        src = ip_transpose(src, true);
    }
    
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    int range = end - start + 1; // +1 because inclusive
    int outputWidth;
    int outputHeight;
    // starts and ends are all inclusive
    int xStart = 0;
    int yStart = 0;
    int xEnd = srcWidth - 1;
    int yEnd = srcHeight - 1;
    outputWidth = srcWidth;
    outputHeight = range;
    yStart = start;
    yEnd = end;
    double energyTable [outputWidth][outputHeight]; // DP energy of pixel
    int parentTable [outputWidth][outputHeight]; // where the previous pixel in the seam came from
    // fill out the energyTable and parentTable
    for (int y = yStart; y <= yEnd; y++) {
        int tableY = y - yStart;
        for (int x = xStart; x <= xEnd; x++) {
            // don't cut a seam through the 'transparent' border
            if (ip_is_border(src, x, y)) {
                energyTable[x][tableY] = INFT * 10;
                parentTable[x][tableY] = x;
                continue;
            }
            if (y == yStart) {
                energyTable[x][tableY] = ip_get_energy(src, x, y);
                parentTable[x][tableY] = -1; // a pixel in the top row has no parent
            } else {
                double energyOfPathFromTop = energyTable[x][tableY - 1] + ip_get_forward_energy(src, x, y, DIRECTION_TOP);
                // the energies of the paths to top left or top right depend on the neighbor existing
                double energyOfPathFromTopLeft = INFT;
                double energyOfPathFromTopRight = INFT;
                if (x > 0) {
                    energyOfPathFromTopLeft = energyTable[x - 1][tableY - 1] + ip_get_forward_energy(src, x, y, DIRECTION_TOP_LEFT);
                }
                if (x + 1 < srcWidth) {
                    energyOfPathFromTopRight = energyTable[x + 1][tableY - 1] + ip_get_forward_energy(src, x, y, DIRECTION_TOP_RIGHT);
                }
                double minParentEnergy = min(energyOfPathFromTop, min(energyOfPathFromTopLeft, energyOfPathFromTopRight));
                energyTable[x][tableY] = ip_get_energy(src, x, y) + minParentEnergy;
                
                if (minParentEnergy == energyOfPathFromTop) {
                    parentTable[x][tableY] = x;
                } else if (minParentEnergy == energyOfPathFromTopLeft) {
                    parentTable[x][tableY] = max(x - 1, 0);
                } else { // energyOfPathFromTopRight
                    parentTable[x][tableY] = x + 1;
                }
            }
        }
    }
    
    int *seam = new int[range];
    // determine the coordinate on the bottom row with the lowest path energy
    double minWeight = INFT;
    double minWeightCoord = 0;
    for (int x = 0; x < srcWidth; x++) {
        if (energyTable[x][range - 1] < minWeight) {
            minWeight = energyTable[x][range - 1];
            minWeightCoord = x;
        }
    }
    seam[range - 1] = minWeightCoord;
    // now that we have the end of the seam determined, backtrack up the parents to find the rest
    for (int i = range - 1 - 1; i >= 0; i--) {
        seam[i] = parentTable[seam[i + 1]][i + 1];
    }
    return seam;
}

/*
 * Computes the lowest-energy seam and returns the coordinates of its pixels.
 */
int* ip_get_seam(Image* src, SeamOrientation orientation)
{
    if (orientation == ORIENTATION_VERTICAL) {
        return ip_get_seam_in_range(src, orientation, 0, src->getHeight() - 1);
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        return ip_get_seam_in_range(src, orientation, 0, src->getWidth() - 1);
    } else {
        // this case should never happen
        return nullptr;
    }
}

/*
 * Creates a new image with the lowest-cost seam shown in red.
 */
Image* ip_show_seam(Image* src, SeamOrientation orientation)
{
    if (orientation == ORIENTATION_HORIZONTAL) {
        src = ip_transpose(src, true);
    }
    int *seam = ip_get_seam(src, ORIENTATION_VERTICAL);
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    Pixel redPixel = Pixel(1, 0, 0);
    
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            if (x == seam[y]) {
                output->setPixel(x, y, redPixel);
            } else {
                output->setPixel(x, y, srcPixel);
            }
        }
    }
    if (orientation == ORIENTATION_HORIZONTAL) {
        return ip_transpose(output, false);
    } else {
        return output;
    }
}

/*
 * Creates an image one pixel less wide or tall with the lowest-cost seam removed.
 */
Image* ip_carve_seam(Image* src, SeamOrientation orientation)
{
    int *seam = ip_get_seam(src, orientation);
    
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    int outputWidth;
    int outputHeight;
    if (orientation == ORIENTATION_VERTICAL) {
        outputWidth = srcWidth - 1;
        outputHeight = srcHeight;
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        outputWidth = srcWidth;
        outputHeight = srcHeight - 1;
    }
    Image* output = new Image(outputWidth, outputHeight);
    
    for (int x = 0; x < outputWidth; x++) {
        for (int y = 0; y < outputHeight; y++) {
            Pixel srcPixel;
            if (orientation == ORIENTATION_VERTICAL) {
                if (x >= seam[y]) {
                    // the current pixel is to the right of the seam (or in it), so sample from the right
                    srcPixel = src->getPixel(x + 1, y);
                } else {
                    // this pixel is to the left of the seam and is thus unaffected by the carving
                    srcPixel = src->getPixel(x, y);
                }
            } else if (orientation == ORIENTATION_HORIZONTAL) {
                if (y >= seam[x]) {
                    // the current pixel is to the bottom of the seam (or in it), so sample from above
                    srcPixel = src->getPixel(x, y + 1);
                } else {
                    // this pixel is above the seam and is thus unaffected by the carving
                    srcPixel = src->getPixel(x, y);
                }
            }
            output->setPixel(x, y, srcPixel);
        }
    }
    return output;
}

/*
 * Creates an image by removing numSeams vertical seams from the given source image.
 */
Image* ip_carve_seams(Image* src, SeamOrientation orientation, int numSeams)
{
    if (orientation == ORIENTATION_VERTICAL) {
        int srcWidth = src->getWidth();
        if (numSeams > srcWidth) {
            cerr << "Cannot carve more seams than the width of the image" << endl;
            return new Image(*src);
        } else {
            Image* resultImage = new Image(*src);
            for (int i = 0; i < numSeams; i++) {
                resultImage = ip_carve_seam(resultImage, orientation);
            }
            return resultImage;
        }
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        int srcHeight = src->getHeight();
        if (numSeams > srcHeight) {
            cerr << "Cannot carve more seams than the height of the image" << endl;
            return new Image(*src);
        } else {
            Image* resultImage = new Image(*src);
            for (int i = 0; i < numSeams; i++) {
                resultImage = ip_carve_seam(resultImage, orientation);
            }
            return resultImage;
        }
    } else {
        // this case should never happen
        return nullptr;
    }
}

/*
 * Returns a new image with the lowest-cost seam duplicated in the given range, shifting over pixels.
 */
Image* ip_insert_seam_in_range(Image* src, SeamOrientation orientation, int start, int end, bool shiftToEnd)
{
    int width = src->getWidth();
    int height = src->getHeight();
    if (orientation == ORIENTATION_VERTICAL) {
        // note: this step wastefully calls the underlying get_seam method twice
        int *seam = ip_get_seam_in_range(src, orientation, start, end);
        Image *outputImage = new Image(width, height);
        for (int x = 0; x < width; x++) {
            // above range
            for (int y = 0; y < start; y++) {
                Pixel srcPixel = src->getPixel(x, y);
                outputImage->setPixel(x, y, srcPixel);
            }
            // below range
            for (int y = end + 1; y < height; y++) {
                Pixel srcPixel = src->getPixel(x, y);
                outputImage->setPixel(x, y, srcPixel);
            }
            // in range
            for (int y = start; y <= end; y++) {
                int seamY = y - start;
                if (x == seam[seamY]) {
                    // average
                    int srcX;
                    if (shiftToEnd) {
                        srcX = x - 1;
                    } else {
                        srcX = x + 1;
                    }
                    if (srcX - 1 >= 0 && srcX + 1 < width) {
                        Pixel leftPixel = src->getPixel(srcX - 1, y);
                        Pixel rightPixel = src->getPixel(srcX + 1, y);
                        Pixel averagePixel = Pixel();
                        averagePixel.setColor(RED, (leftPixel.getColor(RED) + rightPixel.getColor(RED)) / 2);
                        averagePixel.setColor(GREEN, (leftPixel.getColor(GREEN) + rightPixel.getColor(GREEN)) / 2);
                        averagePixel.setColor(BLUE, (leftPixel.getColor(BLUE) + rightPixel.getColor(BLUE)) / 2);
                        outputImage->setPixel(x, y, averagePixel);
                    } else if (srcX - 1 >= 0) {
                        Pixel leftPixel = src->getPixel(srcX - 1, y);
                        outputImage->setPixel(x, y, leftPixel);
                    } else if (srcX + 1 < width) {
                        Pixel rightPixel = src->getPixel(srcX + 1, y);
                        outputImage->setPixel(x, y, rightPixel);
                    } else {
                        Pixel srcPixel = src->getPixel(srcX, y);
                        outputImage->setPixel(x, y, srcPixel);
                    }
                } else {
                    if (x > seam[seamY] && shiftToEnd) {
                        // shift
                        int srcX = x - 1;
                        Pixel srcPixel = src->getPixel(srcX, y);
                        outputImage->setPixel(x, y, srcPixel);
                    } else if (x < seam[seamY] && !shiftToEnd) {
                        // shift
                        int srcX = x + 1;
                        Pixel srcPixel = src->getPixel(srcX, y);
                        outputImage->setPixel(x, y, srcPixel);
                    } else {
                        // do nothing
                        Pixel srcPixel = src->getPixel(x, y);
                        outputImage->setPixel(x, y, srcPixel);
                    }
                }
            }
        }
        return outputImage;
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        // note: this step wastefully calls the underlying get_seam method twice
        int *seam = ip_get_seam_in_range(src, orientation, start, end);
        Image *outputImage = new Image(width, height);
        for (int y = 0; y < height; y++) {
            // left of range
            for (int x = 0; x < start; x++) {
                Pixel srcPixel = src->getPixel(x, y);
                outputImage->setPixel(x, y, srcPixel);
            }
            // right of range
            for (int x = end + 1; x < width; x++) {
                Pixel srcPixel = src->getPixel(x, y);
                outputImage->setPixel(x, y, srcPixel);
            }
            // in range
            for (int x = start; x <= end; x++) {
                int seamX = x - start;
                if (y == seam[seamX]) {
                    // average
                    int srcY;
                    if (shiftToEnd) {
                        srcY = y - 1;
                    } else {
                        srcY = y + 1;
                    }
                    if (srcY - 1 >= 0 && srcY + 1 < height) {
                        Pixel topPixel = src->getPixel(x, srcY - 1);
                        Pixel bottomPixel = src->getPixel(x, srcY + 1);
                        Pixel averagePixel = Pixel();
                        averagePixel.setColor(RED, (topPixel.getColor(RED) + bottomPixel.getColor(RED)) / 2);
                        averagePixel.setColor(GREEN, (topPixel.getColor(GREEN) + bottomPixel.getColor(GREEN)) / 2);
                        averagePixel.setColor(BLUE, (topPixel.getColor(BLUE) + bottomPixel.getColor(BLUE)) / 2);
                        outputImage->setPixel(x, y, averagePixel);
                    } else if (srcY - 1 >= 0) {
                        Pixel topPixel = src->getPixel(x, srcY - 1);
                        outputImage->setPixel(x, y, topPixel);
                    } else if (srcY + 1 < height) {
                        Pixel bottomPixel = src->getPixel(x, srcY + 1);
                        outputImage->setPixel(x, y, bottomPixel);
                    } else {
                        Pixel srcPixel = src->getPixel(x, srcY);
                        outputImage->setPixel(x, y, srcPixel);
                    }
                } else {
                    if (y > seam[seamX] && shiftToEnd) {
                        // shift
                        int srcY = y - 1;
                        Pixel srcPixel = src->getPixel(x, srcY);
                        outputImage->setPixel(x, y, srcPixel);
                    } else if (y < seam[seamX] && !shiftToEnd) {
                        // shift
                        int srcY = y + 1;
                        Pixel srcPixel = src->getPixel(x, srcY);
                        outputImage->setPixel(x, y, srcPixel);
                    } else {
                        // do nothing
                        Pixel srcPixel = src->getPixel(x, y);
                        outputImage->setPixel(x, y, srcPixel);
                    }
                }
            }
        }
        return outputImage;
    } else {
        // this case should never happen
        return nullptr;
    }
}

/*
 * Returns a new image with the lowest-cost seam duplicated.
 */
Image* ip_insert_seam(Image* src, SeamOrientation orientation, bool shiftToEnd)
{
    if (orientation == ORIENTATION_VERTICAL) {
        return ip_insert_seam_in_range(src, orientation, 0, src->getHeight() - 1, shiftToEnd);
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        return ip_insert_seam_in_range(src, orientation, 0, src->getWidth() - 1, shiftToEnd);
    } else {
        // this case should never happen
        return nullptr;
    }
}

/*
 * Inserts numSeams seams into the image by repeatedly inserting one seam. Note that this means all
 * the inserted seams will be the same.
 */
Image* ip_insert_seams(Image* src, SeamOrientation orientation, int numSeams)
{
    Image* resultImage = new Image(*src);
    for (int i = 0; i < numSeams; i++) {
        resultImage = ip_insert_seam(resultImage, orientation, false);
    }
    return resultImage;
}

/*
 * Returns the longest line of 'transparent' pixels along an edge of the image.
 */
vector<int> ip_get_longest_boundary(Image* src, BorderSide& whichSide)
{
    vector<int> boundary;
    size_t maxBoundaryLength = 0;
    int width = src->getWidth();
    int height = src->getHeight();
    
    // top
    bool isCountingBoundary = false;
    vector<int> currentBoundary;
    for (int x = 0; x < width; x++) {
        Pixel leftPixel = src->getPixel(x, 0);
        if (!ip_is_transparent(leftPixel) || (x == width - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_TOP;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(x);
            isCountingBoundary = true;
        }
    }
    
    // bottom
    isCountingBoundary = false;
    currentBoundary.clear();
    for (int x = 0; x < width; x++) {
        Pixel leftPixel = src->getPixel(x, height - 1);
        if (!ip_is_transparent(leftPixel) || (x == width - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_BOTTOM;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(x);
            isCountingBoundary = true;
        }
    }
    
    // left
    isCountingBoundary = false;
    currentBoundary.clear();
    for (int y = 0; y < height; y++) {
        Pixel topPixel = src->getPixel(0, y);
        if (!ip_is_transparent(topPixel) || (y == height - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_LEFT;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(y);
            isCountingBoundary = true;
        }
    }
    
    // right
    isCountingBoundary = false;
    currentBoundary.clear();
    for (int y = 0; y < height; y++) {
        Pixel topPixel = src->getPixel(width - 1, y);
        if (!ip_is_transparent(topPixel) || (y == height - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_RIGHT;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(y);
            isCountingBoundary = true;
        }
    }
    
    return boundary;
}

/*
 * Warps the given image with an irregular 'transparent' border to fill the rectangle.
 * Note that the 'transparent' border must NOT have any boundary as long as its dimension.
 */
Image* ip_local_warp(Image* src)
{
    Image* resultImage = new Image(*src);
    while (true) {
        BorderSide whichSide;
        vector<int> longestBoundary = ip_get_longest_boundary(resultImage, whichSide);
        if (longestBoundary.size() == 0) {
            return resultImage;
        } else {
            SeamOrientation orientation;
            bool shiftToEnd;
            if (whichSide == BORDER_TOP) {
                orientation = ORIENTATION_HORIZONTAL;
                shiftToEnd = false;
            } else if (whichSide == BORDER_BOTTOM) {
                orientation = ORIENTATION_HORIZONTAL;
                shiftToEnd = true;
            } else if (whichSide == BORDER_LEFT) {
                orientation = ORIENTATION_VERTICAL;
                shiftToEnd = false;
            } else if (whichSide == BORDER_RIGHT) {
                orientation = ORIENTATION_VERTICAL;
                shiftToEnd = true;
            }
            int start = longestBoundary[0];
            int end = longestBoundary[longestBoundary.size() - 1];
            resultImage = ip_insert_seam_in_range(resultImage, orientation, start, end, shiftToEnd);
        }
    }
    return resultImage;
}

/*
 * Fits the given src image to its rectangular boundaries, using white as transparent.
 * This is the main method of the project.
 */
Image* ip_rectangle(Image* src)
{
    1; // Currently, this method just replaces all white pixels with blue pixels to show the boundaries.
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image(width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            double red = srcPixel.getColor(RED);
            double green = srcPixel.getColor(GREEN);
            double blue = srcPixel.getColor(BLUE);
            if (red == 1 && green == 1 && blue == 1) {
                Pixel bluePixel = Pixel(0, 0, 1);
                output->setPixel(x, y, bluePixel);
            } else {
                output->setPixel(x, y, srcPixel);
            }
        }
    }
    return output;
}
