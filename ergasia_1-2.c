#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// CAR IMAGE
#define N 278 // height
#define M 420 // width
#define filename "car_420x278_444.yuv"
#define file_yuv "outputCar.yuv"

// CAT
// #define N 332
// #define M 498
// #define filename "cat_498x332_444.yuv"
// #define file_yuv "outputCat.yuv"

// SUNFLOWER
// #define N 200 // height
// #define M 200 // width
// #define filename "sunflower_200x200_444.yuv"
// #define file_yuv "outputFlower.yuv"

/* code for armulator*/
#pragma arm section zidata = "ram"
int current_y[N][M];
int current_u[N][M];
int current_v[N][M];
#pragma arm section

int i, j;

void read_img(void);
void write_img(char *name);

int **gaussianBlur(int **channel, int size, int sigma);
int **sobel(int **channel);
int **nms(int **gradMag, int **gradDir);
int **thresholding(int **channel, int low, int high, int weak);
int **hysteresis(int **channel, int weak);

void copyArray(int **source, int **destination, int rows, int cols);
double **gaussianKernel(int size, int sigma);
void sobelKernelX(double **kernel);
void sobelKernelY(double **kernel);
int **convolution(int **image, int kernelSize, double **kernel);

// If your 2D arrays are stored as a single contiguous block of memory (e.g., a flattened 1D array), you can use memcpy for faster copying:
int **allocate2DArray(int rows, int cols);
void copyFromDynamicToStatic(int **source, int destination[N][M], int rows, int cols);
void copyFromStaticToDynamic(int source[N][M], int **destination, int rows, int cols);

int main()
{
    int **yChannel = allocate2DArray(N, M);
    int **blurredImage = allocate2DArray(N, M);
    int **sobelImage = allocate2DArray(N, M);
    int **gradMag = allocate2DArray(N, M);
    int **gradDir = allocate2DArray(N, M);
    int **nmsImage = allocate2DArray(N, M);
    int **tImage = allocate2DArray(N, M);
    int **hImage = allocate2DArray(N, M);
    int weak = 50;

    read_img();

    copyFromStaticToDynamic(current_y, yChannel, N, M);
    // 1. Gaussian Blur
    blurredImage = gaussianBlur(yChannel, 7, 1);
    copyFromDynamicToStatic(blurredImage, current_y, N, M);
    write_img("BlurredImage.yuv");

    // 2. Sobel Filter
    sobelImage = sobel(blurredImage);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            gradMag[i][j] = sobelImage[i][j];
            gradDir[i][j] = sobelImage[i + N][j];
        }
    }
    copyFromDynamicToStatic(gradMag, current_y, N, M);
    write_img("GradMagImage.yuv");

    // 3. Non-maximum Suppression
    nmsImage = nms(gradMag, gradDir);
    copyFromDynamicToStatic(nmsImage, current_y, N, M);
    write_img("NMSImage.yuv");

    // // 4. Hysterisis Thresholding
    tImage = thresholding(nmsImage, 5, 50, weak);
    copyFromDynamicToStatic(nmsImage, current_y, N, M);
    write_img("ThreshImage.yuv");

    hImage = hysteresis(tImage, weak);
    copyFromDynamicToStatic(hImage, current_y, N, M);

    write_img("FinalImage.yuv");

    return 0;
}

void read_img()
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

void write_img(char *name)
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

    // for (i = 0; i < N; i++)
    // {
    //     for (j = 0; j < M; j++)
    //     {
    //         // fputc(current_u[i][j], frame_yuv);
    //         fputc(0, frame_yuv);
    //     }
    // }

    // for (i = 0; i < N; i++)
    // {
    //     for (j = 0; j < M; j++)
    //     {
    //         // fputc(current_v[i][j], frame_yuv);
    //         fputc(0, frame_yuv);
    //     }
    // }

    fclose(frame_yuv);
}

int **gaussianBlur(int **channel, int size, int sigma)
{
    int kernelSize = size;
    int kernelCenter = kernelSize / 2;
    int kernelSigma = sigma;

    int **copyChannel = allocate2DArray(N, M);

    double **kernel = gaussianKernel(kernelSize, kernelSigma);
    copyArray(channel, copyChannel, N, M);

    // Allocate memory for the blurred channel
    return convolution(channel, kernelSize, kernel);
}

int **sobel(int **channel)
{
    int **output = allocate2DArray(2 * N, M);
    int i, j;
    int **gradX;
    int **gradY;

    // Allocate memory for the Sobel kernel
    double **Gx = malloc(3 * sizeof(double *));
    double **Gy = malloc(3 * sizeof(double *));
    double maxG = 0.0f;
    double G = 0.0f;
    double normCoeff = 0.0f;

    for (i = 0; i < 3; i++)
    {
        Gx[i] = malloc(3 * sizeof(double));
        Gy[i] = malloc(3 * sizeof(double));
    }

    sobelKernelX(Gx);
    sobelKernelY(Gy);
    gradX = convolution(channel, 3, Gx);
    gradY = convolution(channel, 3, Gy);

    // Gradient Magnitude

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            G = sqrt(gradX[i][j] * gradX[i][j] + gradY[i][j] * gradY[i][j]);
            if (G > maxG)
            {
                maxG = G;
            }

            output[i][j] = G;
        }
    }

    normCoeff = 255 / maxG;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            // output[i][j] = (int)fmin(output[i][j] * normCoeff, 255);
            output[i][j] = (output[i][j] * normCoeff > 255) ? 255 : (int)(output[i][j] * normCoeff);
        }
    }

    for (i = N; i < 2 * N; i++)
    {
        for (j = 0; j < M; j++)
        {
            output[i][j] = atan2(gradY[i - N][j], gradX[i - N][j]) * 180 / 3.14159265359;
        }
    }

    for (i = 0; i < 3; i++)
    {
        free(Gx[i]);
        free(Gy[i]);
    }
    free(Gx);
    free(Gy);

    return output;
}

int **nms(int **gradMag, int **gradDir)
{
    int **output = allocate2DArray(N, M);
    int dir, beforePixel, afterPixel;
    int PI = 180;
    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            dir = gradDir[i][j];

            if ((0 <= dir < PI / 8) || (15 * PI / 8 <= dir <= 2 * PI))
            {
                beforePixel = gradMag[i][j - 1];
                afterPixel = gradMag[i][j + 1];
            }
            else if ((PI / 8 <= dir < 3 * PI / 8) || (9 * PI / 8 <= dir <= 11 * PI / 8))
            {
                beforePixel = gradMag[i + 1][j - 1];
                afterPixel = gradMag[i - 1][j + 1];
            }
            else if ((3 * PI / 8 <= dir < 5 * PI / 8) || (11 * PI / 8 <= dir <= 13 * PI / 8))
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
        }
    }

    return output;
}

int **thresholding(int **channel, int low, int high, int weak)
{
    int **output = allocate2DArray(N, M);

    int strong = 255;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
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
    int **output = allocate2DArray(N, M);
    int **topToBottom = allocate2DArray(N, M);
    int **bottomToTop = allocate2DArray(N, M);
    int **rightToLeft = allocate2DArray(N, M);
    int **leftToRight = allocate2DArray(N, M);
    // int topToBottom[N][M], bottomToTop[N][M], rightToLeft[N][M], leftToRight[N][M];
    int fp = -1;

    copyArray(channel, topToBottom, N, M);
    copyArray(channel, bottomToTop, N, M);
    copyArray(channel, rightToLeft, N, M);
    copyArray(channel, leftToRight, N, M);

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

    return output;
}

double **gaussianKernel(int size, int sigma)
{
    double **kernel = malloc(sizeof(double *) * size);
    int i, j;
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

    for (i = 0; i < size; i++) // checking purpose
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

void copyArray(int **source, int **destination, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            destination[i][j] = source[i][j];
        }
    }
}

int **convolution(int **image, int kernelSize, double **kernel)
{
    int i, j, ki, kj;
    int kernelCenter = kernelSize / 2;
    int sum = 0;

    int **output = malloc(N * sizeof(int *));

    for (i = 0; i < N; i++)
    {
        output[i] = malloc(M * sizeof(int));
        for (j = 0; j < M; j++)
        {
            sum = 0;
            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
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

    return output;
}

int **allocate2DArray(int rows, int cols)
{
    // Allocate memory for the row pointers
    int **array = (int **)malloc(rows * sizeof(int *));
    if (array == NULL)
    {
        perror("Failed to allocate memory for rows");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for each row
    for (int i = 0; i < rows; i++)
    {
        array[i] = (int *)malloc(cols * sizeof(int));
        if (array[i] == NULL)
        {
            perror("Failed to allocate memory for columns");
            // Free already allocated rows
            for (int j = 0; j < i; j++)
            {
                free(array[j]);
            }
            free(array);
            exit(EXIT_FAILURE);
        }
    }

    return array;
}

void copyFromDynamicToStatic(int **source, int destination[N][M], int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            destination[i][j] = source[i][j];
        }
    }
}

void copyFromStaticToDynamic(int source[N][M], int **destination, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            destination[i][j] = source[i][j];
        }
    }
}