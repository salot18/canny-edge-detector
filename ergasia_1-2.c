#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* CAR IMAGE */
// #define N 278 // Image height (rows)
// #define M 420 // Image width (cols)
// #define filename "car_420x278_444.yuv"

/* CAT IMAGES */
// #define N 332
// #define M 498
// #define filename "cat_498x332_444.yuv"

// SUNFLOWER
#define N 200 // height
#define M 200 // width
#define filename "sunflower_200x200_444.yuv"

#define _PI 3.14159265359

/* code for armulator*/
#pragma arm section zidata = "ram"
int current_y[N][M];
int current_u[N][M];
int current_v[N][M];
#pragma arm section

int i, j, ii, jj;
int TILE_SIZE = 100;

/* FUNCTIONS */
void readImage(void);
void writeImage(char *name);

/* Canny Algorithm */
void gaussianBlur(int **channel, int size, int sigma, int **output);
int **sobel(int **channel);
int **nms(int **gradMag, int **gradDir);
int **thresholding(int **channel, int low, int high, int weak);
int **hysteresis(int **channel, int weak);

double **gaussianKernel(int size, int sigma);
void sobelKernelX(double **kernel);
void sobelKernelY(double **kernel);
void convolution(int **image, int kernelSize, double **kernel, int **output);

/* Array Memory Management */
int **allocate2DIntArray(int rows, int cols);
double **allocate2DDoubleArray(int rows, int cols);
void freeInt2DArray(int **array, int rows, int cols);
void freeDouble2DArray(double **array, int rows, int cols);

void copyFromDynamicToDynamic(int **source, int **destination, int rows, int cols);
void copyFromDynamicToStatic(int **source, int destination[N][M], int rows, int cols);
void copyFromStaticToDynamic(int source[N][M], int **destination, int rows, int cols);

void thresholdCheck(int *channel, int *output, int low, int high, int weak, int strong);

int main()
{
    int **yChannel = allocate2DIntArray(N, M);
    int **blurredImage = allocate2DIntArray(N, M);
    int **sobelImage;
    int **gradMag = allocate2DIntArray(N, M);
    int **gradDir = allocate2DIntArray(N, M);
    int **nmsImage;
    int **tImage;
    int **hImage;
    int weak = 50;

    readImage();

    copyFromStaticToDynamic(current_y, yChannel, N, M);

    /* 1. GAUSSIAN BLUR */
    gaussianBlur(yChannel, 7, 1, blurredImage); // (image, kernelSize, kernelSigma, output)
    // copyFromDynamicToStatic(blurredImage, current_y, N, M);  // Get the blurred image
    // writeImage("BlurredImage.yuv");                            // and save it

    /* 2. SOBEL MASK */
    sobelImage = sobel(blurredImage);
    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj++)
                {
                    gradMag[ii][jj] = sobelImage[ii][jj];
                    gradDir[ii][jj] = sobelImage[ii + N][jj];
                }
            }
        }
    }
    // copyFromDynamicToStatic(gradMag, current_y, N, M);      // Get the gradient magnitude image
    // writeImage("GradMagImage.yuv");                           // and save it

    /* 3. NON-MAXIMUM SUPPRESSION */
    nmsImage = nms(gradMag, gradDir);
    // copyFromDynamicToStatic(nmsImage, current_y, N, M);     // Get the image afte the nms
    // writeImage("NMSImage.yuv");                               // and save it

    /* 4. HYSTERESIS THRESHOLDING */
    tImage = thresholding(nmsImage, 5, 20, weak); // (image, low, hight, weak)
    // copyFromDynamicToStatic(tImage, current_y, N, M);       // Get the image after thresholding
    // writeImage("ThreshImage.yuv");                            // and save it
    hImage = hysteresis(tImage, weak);
    copyFromDynamicToStatic(hImage, current_y, N, M); // Get the image after hysteresis (final image)
    writeImage("FinalImage.yuv");                     // and save it

    /* Free memory for all dynamically allocated 2D arrays */
    freeInt2DArray(yChannel, N, M);
    freeInt2DArray(blurredImage, N, M);
    freeInt2DArray(sobelImage, 2 * N, M);
    freeInt2DArray(gradMag, N, M);
    freeInt2DArray(gradDir, N, M);
    freeInt2DArray(nmsImage, N, M);
    freeInt2DArray(tImage, N, M);
    freeInt2DArray(hImage, N, M);

    return 0;
}

void readImage(void)
{
    FILE *frame_c;
    if ((frame_c = fopen(filename, "rb")) == NULL)
    {
        printf("current frame doesn't exist\n");
        exit(-1);
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            current_y[i][j] = fgetc(frame_c);
            current_y[i][j + 1] = fgetc(frame_c);
            current_y[i][j + 2] = fgetc(frame_c);
            current_y[i][j + 3] = fgetc(frame_c);
        }
    }

    fclose(frame_c);
}

void writeImage(char *name)
{
    FILE *frame_yuv;
    frame_yuv = fopen(name, "wb");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            fputc(current_y[i][j], frame_yuv);
            fputc(current_y[i][j + 1], frame_yuv);
            fputc(current_y[i][j + 2], frame_yuv);
            fputc(current_y[i][j + 3], frame_yuv);
        }
    }

    fclose(frame_yuv);
}

void gaussianBlur(int **channel, int kernelSize, int kernelSigma, int **output)
{
    // Generate the Gaussian kernel for size and sigma
    double **kernel = gaussianKernel(kernelSize, kernelSigma);

    convolution(channel, kernelSize, kernel, output);

    freeDouble2DArray(kernel, kernelSize, kernelSize);
}

int **sobel(int **channel)
{
    double maxG = 0.0f; // Maximum gradient magnitude value, used for normalization
    double G = 0.0f;    // Gradient Magnitude
    double normCoeff = 0.0f;

    // Allocate memory for the output, the Sobel kernels and the directional gradients
    int **output = allocate2DIntArray(2 * N, M);
    int **gradX = allocate2DIntArray(N, M);
    int **gradY = allocate2DIntArray(N, M);
    double **Gx = allocate2DDoubleArray(3, 3);
    double **Gy = allocate2DDoubleArray(3, 3);

    sobelKernelX(Gx);
    sobelKernelY(Gy);

    convolution(channel, 3, Gx, gradX);
    convolution(channel, 3, Gy, gradY);

    // Gradient Magnitude
    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {
                    output[ii][jj + 0] = sqrt(gradX[ii][jj + 0] * gradX[ii][jj + 0] + gradY[ii][jj + 0] * gradY[ii][jj + 0]);
                    output[ii][jj + 1] = sqrt(gradX[ii][jj + 1] * gradX[ii][jj + 1] + gradY[ii][jj + 1] * gradY[ii][jj + 1]);
                    output[ii][jj + 2] = sqrt(gradX[ii][jj + 2] * gradX[ii][jj + 2] + gradY[ii][jj + 2] * gradY[ii][jj + 2]);
                    output[ii][jj + 3] = sqrt(gradX[ii][jj + 3] * gradX[ii][jj + 3] + gradY[ii][jj + 3] * gradY[ii][jj + 3]);
                }
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            if (output[i][j] > maxG)
            {
                maxG = output[i][j];
            }
        }
    }

    // Normalize gradient magnitude (0 - 255)
    normCoeff = 255 / maxG;
    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {
                    output[ii][jj + 0] = (output[ii][jj + 0] * normCoeff > 255) ? 255 : (int)(output[ii][jj + 0] * normCoeff);
                    output[ii][jj + 1] = (output[ii][jj + 1] * normCoeff > 255) ? 255 : (int)(output[ii][jj + 1] * normCoeff);
                    output[ii][jj + 2] = (output[ii][jj + 2] * normCoeff > 255) ? 255 : (int)(output[ii][jj + 2] * normCoeff);
                    output[ii][jj + 3] = (output[ii][jj + 3] * normCoeff > 255) ? 255 : (int)(output[ii][jj + 3] * normCoeff);
                }
            }
        }
    }

    // Gradient Direction
    for (i = N; i < 2 * N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {

                    output[i][j + 0] = atan2(gradY[i - N][j + 0], gradX[i - N][j + 0]) * 180 / _PI;
                    output[i][j + 1] = atan2(gradY[i - N][j + 1], gradX[i - N][j + 1]) * 180 / _PI;
                    output[i][j + 2] = atan2(gradY[i - N][j + 2], gradX[i - N][j + 2]) * 180 / _PI;
                    output[i][j + 3] = atan2(gradY[i - N][j + 3], gradX[i - N][j + 3]) * 180 / _PI;
                }
            }
        }
    }

    freeDouble2DArray(Gx, 3, 3);
    freeDouble2DArray(Gy, 3, 3);
    freeInt2DArray(gradX, N, M);
    freeInt2DArray(gradY, N, M);

    return output;
}

// TODO: Loop carried dependency
int **nms(int **gradMag, int **gradDir)
{
    int **output = allocate2DIntArray(N, M);
    int dir, beforePixel, afterPixel;
    int PI = 180;

    // Ignore the border pixels
    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            dir = gradDir[i][j];

            // ################## REVIEW THE CONDITIONS #########################

            if (((0 <= dir) || (dir < PI / 8)) || ((15 * PI / 8 <= dir) || (dir <= 2 * PI)))
            {
                beforePixel = gradMag[i][j - 1];
                afterPixel = gradMag[i][j + 1];
            }
            else if (((PI / 8 <= dir) || (dir < 3 * PI / 8)) || ((9 * PI / 8 <= dir) || (dir <= 11 * PI / 8)))
            {
                beforePixel = gradMag[i + 1][j - 1];
                afterPixel = gradMag[i - 1][j + 1];
            }
            else if (((3 * PI / 8 <= dir) || (dir < 5 * PI / 8)) || ((11 * PI / 8 <= dir) || (dir <= 13 * PI / 8)))
            {
                beforePixel = gradMag[i - 1][j];
                afterPixel = gradMag[i + 1][j];
            }
            else
            {
                beforePixel = gradMag[i - 1][j - 1];
                afterPixel = gradMag[i + 1][j + 1];
            }

            if (gradMag[i][j] >= beforePixel && gradMag[i][j] >= afterPixel)
            {
                output[i][j] = gradMag[i][j];
            }
            else
            {
                output[i][j] = 0;
            }
        }
    }

    return output;
}

int **thresholding(int **channel, int low, int high, int weak)
{
    int **output = allocate2DIntArray(N, M);
    int strong = 255;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
        }
    }

    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {
                    thresholdCheck(&channel[ii][jj + 0], &output[ii][jj + 0], low, high, weak, strong);
                    thresholdCheck(&channel[ii][jj + 1], &output[ii][jj + 1], low, high, weak, strong);
                    thresholdCheck(&channel[ii][jj + 2], &output[ii][jj + 2], low, high, weak, strong);
                    thresholdCheck(&channel[ii][jj + 3], &output[ii][jj + 3], low, high, weak, strong);
                }
            }
        }
    }

    return output;
}

int **hysteresis(int **channel, int weak)
{
    int **output = allocate2DIntArray(N, M);
    int **topToBottom = allocate2DIntArray(N, M);
    int **bottomToTop = allocate2DIntArray(N, M);
    int **rightToLeft = allocate2DIntArray(N, M);
    int **leftToRight = allocate2DIntArray(N, M);
    // int topToBottom[N][M], bottomToTop[N][M], rightToLeft[N][M], leftToRight[N][M];
    int fp = -1;

    copyFromDynamicToDynamic(channel, topToBottom, N, M);
    copyFromDynamicToDynamic(channel, bottomToTop, N, M);
    copyFromDynamicToDynamic(channel, rightToLeft, N, M);
    copyFromDynamicToDynamic(channel, leftToRight, N, M);

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (topToBottom[i][j] == weak)
            {
                if (topToBottom[i][j + 1] == 255 || topToBottom[i][j - 1] == 255 ||
                    topToBottom[i - 1][j] == 255 || topToBottom[i + 1][j] == 255 ||
                    topToBottom[i - 1][j - 1] == 255 || topToBottom[i + 1][j - 1] == 255 ||
                    topToBottom[i - 1][j + 1] == 255 || topToBottom[i + 1][j + 1] == 255)
                {
                    topToBottom[i][j] = 255;
                }
                else
                {
                    topToBottom[i][j] = 0;
                }
            }
        }
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (bottomToTop[i][j] == weak)
            {
                if (bottomToTop[i][j + 1] == 255 || bottomToTop[i][j - 1] == 255 ||
                    bottomToTop[i - 1][j] == 255 || bottomToTop[i + 1][j] == 255 ||
                    bottomToTop[i - 1][j - 1] == 255 || bottomToTop[i + 1][j - 1] == 255 ||
                    bottomToTop[i - 1][j + 1] == 255 || bottomToTop[i + 1][j + 1] == 255)
                {

                    bottomToTop[i][j] = 255;
                }
                else
                {
                    bottomToTop[i][j] = 0;
                }
            }
        }
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (rightToLeft[i][j] == weak)
            {
                if (rightToLeft[i][j + 1] == 255 || rightToLeft[i][j - 1] == 255 ||
                    rightToLeft[i - 1][j] == 255 || rightToLeft[i + 1][j] == 255 ||
                    rightToLeft[i - 1][j - 1] == 255 || rightToLeft[i + 1][j - 1] == 255 ||
                    rightToLeft[i - 1][j + 1] == 255 || rightToLeft[i + 1][j + 1] == 255)
                {

                    rightToLeft[i][j] = 255;
                }
                else
                {

                    rightToLeft[i][j] = 0;
                }
            }
        }
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (leftToRight[i][j] == weak)
            {
                if (leftToRight[i][j + 1] == 255 || leftToRight[i][j - 1] == 255 ||
                    leftToRight[i - 1][j] == 255 || leftToRight[i + 1][j] == 255 ||
                    leftToRight[i - 1][j - 1] == 255 || leftToRight[i + 1][j - 1] == 255 ||
                    leftToRight[i - 1][j + 1] == 255 || leftToRight[i + 1][j + 1] == 255)
                {

                    leftToRight[i][j] = 255;
                }
                else
                {

                    leftToRight[i][j] = 0;
                }
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            fp = topToBottom[i][j] + bottomToTop[i][j] + rightToLeft[i][j] + leftToRight[i][j];
            if (fp > 255)
            {
                fp = 255;
            }

            output[i][j] = fp;
        }
    }

    freeInt2DArray(topToBottom, N, M);
    freeInt2DArray(bottomToTop, N, M);
    freeInt2DArray(rightToLeft, N, M);
    freeInt2DArray(leftToRight, N, M);

    return output;
}

// TODO: Calculate only 1/4 of the kernel, due to symmetry
double **gaussianKernel(int size, int sigma)
{
    double **kernel = malloc(sizeof(double *) * size);
    int x, y;
    int center = size / 2;
    double res = 0.0f;
    double sum = 0.0f;

    for (i = 0; i < size; i++) // checking purpose
    {
        kernel[i] = malloc(sizeof(double) * size);
        for (j = 0; j < size; j++)
        {
            x = j - center;
            y = i - center;
            res = exp(-((double)(x * x + y * y) / (2 * sigma * sigma)));
            sum += res;
            kernel[i][j] = res;
        }
    }

    // Normalize the kernel
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            kernel[i][j] = kernel[i][j] / sum;
        }
    }

    return kernel;
}

void sobelKernelX(double **kernel)
{
    kernel[0][0] = -1;
    kernel[0][1] = 0;
    kernel[0][2] = 1;
    kernel[1][0] = -2;
    kernel[1][1] = 0;
    kernel[1][2] = 2;
    kernel[2][0] = -1;
    kernel[2][1] = 0;
    kernel[2][2] = 1;
}

void sobelKernelY(double **kernel)
{
    kernel[0][0] = -1;
    kernel[0][1] = -2;
    kernel[0][2] = -1;
    kernel[1][0] = 0;
    kernel[1][1] = 0;
    kernel[1][2] = 0;
    kernel[2][0] = 1;
    kernel[2][1] = 2;
    kernel[2][2] = 1;
}

void convolution(int **image, int kernelSize, double **kernel, int **output)
{
    int ki, kj;
    int kernelCenter = kernelSize / 2;
    int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            sum0 = sum1 = sum2 = sum3 = 0;
            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
                {
                    if (i + ki >= 0 && i + ki < N && j + kj >= 0 && j + kj < M)
                    {
                        sum0 += image[i + ki][j + kj] * kernel[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 1 + kj >= 0 && j + 1 + kj < M)
                    {
                        sum1 += image[i + ki][j + 1 + kj] * kernel[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 2 + kj >= 0 && j + 2 + kj < M)
                    {
                        sum2 += image[i + ki][j + 2 + kj] * kernel[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 3 + kj >= 0 && j + 3 + kj < M)
                    {
                        sum3 += image[i + ki][j + 3 + kj] * kernel[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            output[i][j] = sum0;
            output[i][j + 1] = sum1;
            output[i][j + 2] = sum2;
            output[i][j + 3] = sum3;
        }
    }

    // return output;
}

int **allocate2DIntArray(int rows, int cols)
{
    // Allocate memory for the row pointers
    int **array = (int **)malloc(rows * sizeof(int *));
    if (array == NULL)
    {
        perror("Failed to allocate memory for rows");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each row
    for (i = 0; i < rows; i++)
    {
        array[i] = (int *)malloc(cols * sizeof(int));
        if (array[i] == NULL)
        {
            perror("Failed to allocate memory for columns");
            // Free already allocated rows
            for (j = 0; j < i; j++)
            {
                free(array[j]);
            }
            free(array);
            exit(EXIT_FAILURE);
        }
    }

    return array;
}

double **allocate2DDoubleArray(int rows, int cols)
{
    // Allocate memory for the row pointers
    double **array = (double **)malloc(rows * sizeof(double *));
    if (array == NULL)
    {
        perror("Failed to allocate memory for rows");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each row
    for (i = 0; i < rows; i++)
    {
        array[i] = (double *)malloc(cols * sizeof(double));
        if (array[i] == NULL)
        {
            perror("Failed to allocate memory for columns");
            // Free already allocated rows
            for (j = 0; j < i; j++)
            {
                free(array[j]);
            }
            free(array);
            exit(EXIT_FAILURE);
        }
    }

    return array;
}

void freeInt2DArray(int **array, int rows, int cols)
{
    for (i = 0; i < rows; i++)
    {
        free(array[i]);
        // free(array[i + 1]);
        // free(array[i + 2]);
        // free(array[i + 3]);
    }
    free(array);
}

void freeDouble2DArray(double **array, int rows, int cols)
{
    for (i = 0; i < rows; i++)
    {
        free(array[i]);
        // free(array[i + 1]);
        // free(array[i + 2]);
        // free(array[i + 3]);
    }
    free(array);
}

void copyFromDynamicToDynamic(int **source, int **destination, int rows, int cols)
{
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j += 4)
        {
            destination[i][j] = source[i][j];
            destination[i][j + 1] = source[i][j + 1];
            destination[i][j + 2] = source[i][j + 2];
            destination[i][j + 3] = source[i][j + 3];
        }
    }
}

void copyFromDynamicToStatic(int **source, int destination[N][M], int rows, int cols)
{
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j += 4)
        {
            destination[i][j] = source[i][j];
            destination[i][j + 1] = source[i][j + 1];
            destination[i][j + 2] = source[i][j + 2];
            destination[i][j + 3] = source[i][j + 3];
        }
    }
}

void copyFromStaticToDynamic(int source[N][M], int **destination, int rows, int cols)
{
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j += 4)
        {
            destination[i][j] = source[i][j];
            destination[i][j + 1] = source[i][j + 1];
            destination[i][j + 2] = source[i][j + 2];
            destination[i][j + 3] = source[i][j + 3];
        }
    }
}

void thresholdCheck(int *channel, int *output, int low, int high, int weak, int strong)
{
    if (*channel >= high)
    {
        *output = strong;
    }
    else if (*channel >= low)
    {
        *output = weak;
    }
    else
    {
        *output = 0;
    }
}