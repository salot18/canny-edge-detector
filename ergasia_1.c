#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// CAR IMAGE
// #define N 278 // height
// #define M 420 // width
// #define filename "./444/car_420x278_444.yuv"
// #define file_yuv "outputCar.yuv"

// CAT
#define N 332
#define M 498
#define filename "./444/cat_498x332_444.yuv"
#define file_yuv "outputCat.yuv"

// SUNFLOWER
// #define N 200 // height
// #define M 200 // width
// #define filename "./444/sunflower_200x200_444.yuv"
// #define file_yuv "outputFlower.yuv"

/* code for armulator*/
#pragma arm section zidata = "ram"
int current_y[N][M];
int current_u[N][M];
int current_v[N][M];
#pragma arm section

int i, j;

void read_img(void);
void write_img(void);

void gaussianBlur(int channel[N][M], int size, int sigma, int output[N][M]);
void sobel(int channel[N][M], int gradMag[N][M], double gradDir[N][M]);
void nms(int gradMag[N][M], double gradDir[N][M], int output[N][M]);
void thresholding(int channel[N][M], int low, int high, int weak, int output[N][M]);
void hysteresis(int channel[N][M], int weak, int output[N][M]);

void copyArray(int source[N][M], int destination[N][M]);
double **gaussianKernel(int size, int sigma);
void sobelKernelX(double **kernel);
void sobelKernelY(double **kernel);
void convolution(int image[N][M], int kernelSize, double **kernel, int output[N][M]);

int main()
{
    int blurredImage[N][M];
    int gradMag[N][M];
    double gradDir[N][M];
    int nmsImage[N][M];
    int weak = 50;
    int tImage[N][M];
    int hImage[N][M];

    read_img();

    // 1. Gaussian Blur
    gaussianBlur(current_y, 7, 1, blurredImage);
    // copyArray(blurredImage, current_y);

    // 2. Sobel Filter
    sobel(blurredImage, gradMag, gradDir);
    // copyArray(gradMag, current_y);

    // 3. Non-maximum Suppression
    nms(gradMag, gradDir, nmsImage);
    // copyArray(nmsImage, current_y);

    // 4. Hysterisis Thresholding
    thresholding(nmsImage, 5, 50, weak, tImage);
    hysteresis(tImage, weak, hImage);
    copyArray(hImage, current_y);

    write_img();

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

void write_img()
{
    FILE *frame_yuv;
    frame_yuv = fopen(file_yuv, "wb");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            fputc(current_y[i][j], frame_yuv);
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            fputc(current_u[i][j], frame_yuv);
            // fputc(0, frame_yuv);
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            fputc(current_v[i][j], frame_yuv);
            // fputc(0, frame_yuv);
        }
    }

    fclose(frame_yuv);
}

void gaussianBlur(int channel[N][M], int size, int sigma, int output[N][M])
{
    int kernelSize = size;
    int kernelCenter = kernelSize / 2;
    int kernelSigma = sigma;

    int copyChannel[N][M];

    double **kernel = gaussianKernel(kernelSize, kernelSigma);
    copyArray(channel, copyChannel);
    // Allocate memory for the blurred channel
    convolution(channel, kernelSize, kernel, output);
}

void sobel(int channel[N][M], int gradMag[N][M], double gradDir[N][M])
{
    int i, j;
    int gradX[N][M];
    int gradY[N][M];

    // Allocate memory for the Sobel kernel
    double **Gx = malloc(3 * sizeof(double *));
    double **Gy = malloc(3 * sizeof(double *));
    double maxG = 0.0f;
    double G = 0.0f;
    double normCoeff = 255.0 / maxG;

    for (i = 0; i < 3; i++)
    {
        Gx[i] = malloc(3 * sizeof(double));
        Gy[i] = malloc(3 * sizeof(double));
    }

    sobelKernelX(Gx);
    sobelKernelY(Gy);
    convolution(channel, 3, Gx, gradX);
    convolution(channel, 3, Gy, gradY);

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

            gradMag[i][j] = G;
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            gradMag[i][j] = (int)fmin(gradMag[i][j] * normCoeff, 255);
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            gradDir[i][j] = atan2(gradY[i][j], gradX[i][j]) * 180 / 3.14159265359;
        }
    }

    for (i = 0; i < 3; i++)
    {
        free(Gx[i]);
        free(Gy[i]);
    }
    free(Gx);
    free(Gy);
}

void nms(int gradMag[N][M], double gradDir[N][M], int output[N][M])
{
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
}

void thresholding(int channel[N][M], int low, int high, int weak, int output[N][M])
{
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
}

void hysteresis(int channel[N][M], int weak, int output[N][M])
{

    int topToBottom[N][M], bottomToTop[N][M], rightToLeft[N][M], leftToRight[N][M];
    int fp = -1;
    copyArray(channel, topToBottom);
    copyArray(channel, bottomToTop);
    copyArray(channel, rightToLeft);
    copyArray(channel, leftToRight);

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

void copyArray(int source[N][M], int destination[N][M])
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            destination[i][j] = source[i][j];
        }
    }
}

void convolution(int image[N][M], int kernelSize, double **kernel, int output[N][M])
{
    int i, j, ki, kj;
    int kernelCenter = kernelSize / 2;
    int sum = 0;

    for (i = 0; i < N; i++)
    {
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
}