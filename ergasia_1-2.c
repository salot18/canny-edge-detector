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

/* code for armulator*/
#pragma arm section zidata = "ram"
int current_y[N][M];
int current_u[N][M];
int current_v[N][M];
#pragma arm section

int i, j;

/* FUNCTIONS */
void readImage();
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
    for (j = 0; j < M; j++)
    {
        for (i = 0; i < N; i++)
        {
            gradMag[i][j] = sobelImage[i][j];
            gradDir[i][j] = sobelImage[i + N][j];
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

void readImage()
{
    FILE *frame_c;
    if ((frame_c = fopen(filename, "rb")) == NULL)
    {
        printf("current frame doesn't exist\n");
        exit(-1);
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            current_y[i][j] = fgetc(frame_c);
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            current_u[i][j] = fgetc(frame_c);
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            current_v[i][j] = fgetc(frame_c);
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
        for (j = 0; j < M; j++)
        {
            fputc(current_y[i][j], frame_yuv);
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
    for (j = 0; j < M; j++)
    {
        for (i = 0; i < N; i++)
        {
            G = sqrt(gradX[i][j] * gradX[i][j] + gradY[i][j] * gradY[i][j]);
            if (G > maxG)
            {
                maxG = G;
            }

            output[i][j] = G;
        }
    }

    // Normalize gradient magnitude (0 - 255)
    normCoeff = 255 / maxG;
    for (j = 0; j < M; j++)
    {
        for (i = 0; i < N; i++)
        {
            output[i][j] = (output[i][j] * normCoeff > 255) ? 255 : (int)(output[i][j] * normCoeff);
        }
    }

    // Gradient Direction
    for (j = 0; j < M; j++)
    {
        for (i = N; i < 2 * N; i++)
        {
            output[i][j] = atan2(gradY[i - N][j], gradX[i - N][j]) * 180 / 3.14159265359;
        }
    }

    freeDouble2DArray(Gx, 3, 3);
    freeDouble2DArray(Gy, 3, 3);
    freeInt2DArray(gradX, N, M);
    freeInt2DArray(gradY, N, M);

    return output;
}

int **nms(int **gradMag, int **gradDir)
{
    int **output = allocate2DIntArray(N, M);
    int dir, beforePixel, afterPixel;
    int PI = 180;

    // Ignore the border pixels
    for (j = 1; j < M - 1; j++)
    {
        for (i = 1; i < N - 1; i++)
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

    for (j = 0; j < M; j++)
    {
        for (i = 0; i < N; i++)
        {
            if (channel[i][j] >= high)
            {
                output[i][j] = strong;
            }
            else if (channel[i][j] >= low)
            {
                output[i][j] = weak;
            }
            else
            {
                output[i][j] = 0;
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

    for (j = 1; j < M - 1; j++)
    {
        for (i = 1; i < N - 1; i++)
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

    for (j = 1; j < M - 1; j++)
    {
        for (i = 1; i < N - 1; i++)
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

    for (j = 1; j < M - 1; j++)
    {
        for (i = 1; i < N - 1; i++)
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

    for (j = 1; j < M - 1; j++)
    {
        for (i = 1; i < N - 1; i++)
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

    for (j = 0; j < M; j++)
    {
        for (i = 0; i < N; i++)
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
    for (j = 0; j < size; j++)
    {
        for (i = 0; i < size; i++)
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
    int i, j, ki, kj;
    int kernelCenter = kernelSize / 2;
    int sum = 0;

    // int **output = malloc(N * sizeof(int *));

    for (j = 0; j < M; j++)
    {
        for (i = 0; i < N; i++)
        // output[i] = malloc(M * sizeof(int));
        {
            sum = 0;
            for (kj = -kernelCenter; kj <= kernelCenter; kj++)
            {
                for (ki = -kernelCenter; ki <= kernelCenter; ki++)
                {
                    if (i + ki >= 0 && i + ki < N && j + kj >= 0 && j + kj < M)
                    {
                        sum += image[i + ki][j + kj] * kernel[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            output[i][j] = sum;
        }
    }

    // return output;
}

int **allocate2DIntArray(int rows, int cols)
{
    int i, j;
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
    int i, j;
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
    int i, j;
    for (i = 0; i < rows; i++)
    {
        free(array[i]);
    }
    free(array);
}

void freeDouble2DArray(double **array, int rows, int cols)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        free(array[i]);
    }
    free(array);
}

void copyFromDynamicToDynamic(int **source, int **destination, int rows, int cols)
{
    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            destination[i][j] = source[i][j];
        }
    }
}

void copyFromDynamicToStatic(int **source, int destination[N][M], int rows, int cols)
{
    int i, j;
    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            destination[i][j] = source[i][j];
        }
    }
}

void copyFromStaticToDynamic(int source[N][M], int **destination, int rows, int cols)
{
    int i, j;
    for (j = 0; j < cols; j++)
    {
        for (i = 0; i < rows; i++)
        {
            destination[i][j] = source[i][j];
        }
    }
}