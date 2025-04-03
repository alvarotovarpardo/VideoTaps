#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <opencv2/opencv.hpp>
#include <cstring>
#include <fstream>
#include <tuple>
using namespace cv;
using namespace std;

tuple<int, int, char, int, int, char> readTap(const string& tapType) {
    size_t xPos = tapType.find('X');
    if (xPos == string::npos || xPos == 0) {
        cerr << "Invalid format" << endl;
        return {-1, -1, '\0', -1, -1, '\0'};
    }

    int R = stoi(tapType.substr(0, xPos));
    int T = 1;
    char L = '\0';

    size_t posL = xPos + 1;
    while (posL < tapType.length() && isdigit(tapType[posL])) {
        posL++;
    }

    if (posL > xPos + 1) {
        T = stoi(tapType.substr(xPos + 1, posL - (xPos + 1)));
    }

    if (posL < tapType.length() && isalpha(tapType[posL])) {
        L = tapType[posL];
        posL++;
    }

    // Si hay un segundo tap (formato '-1Y4' o similar)
    int R2 = 1, T2 = 1;
    char L2 = '\0';

    if (posL < tapType.length() && tapType[posL] == '-') {
        size_t yPos = tapType.find('Y', posL);
        if (yPos != string::npos) {
            R2 = stoi(tapType.substr(posL + 1, yPos - (posL + 1)));
            size_t posL2 = yPos + 1;
            while (posL2 < tapType.length() && isdigit(tapType[posL2])) {
                posL2++;
            }

            if (posL2 > yPos + 1) {
                T2 = stoi(tapType.substr(yPos + 1, posL2 - (yPos + 1)));
            }

            if (posL2 < tapType.length() && isalpha(tapType[posL2])) {
                L2 = tapType[posL2];
            }
        }
    }

    return {R, T, L, R2, T2, L2};
}


void applyTap(const uchar* input, uchar* output, int rows, int cols, const string& tapType) {
    auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);

    int regionWidth = cols / R;
    int tapsNumber = regionWidth / T; 
    int regionHeight = rows / Ry;
    int tapsNumberY = regionHeight / Ty;
    
    
    for (int i = 0; i < regionHeight; i++){
        for (int r = 0; r < R; r++){
            for (int j = 0; j < tapsNumber; j++){
                for (int t = 0; t < T; t++){
                    int srcIndex = (i * cols) + (r * regionWidth + j * T + t);
                    int dstIndex = (i * cols) + (T * (j * R + r) + t);
                    output[dstIndex] = input[srcIndex];
                }
            }
        }
    } 
    /*
    for (int ry = 0; ry < Ry; ry++){
        for (int i = 0; i < tapsNumberY; i++){
            for (int ty = 0; ty < Ty; ty++){
                for (int r = 0; r < R; r++){
                    for (int j = 0; j < tapsNumber; j++){
                        for (int t = 0; t < T; t++){
                            int srcIndex = ((ry * regionHeight + i * tapsNumberY + ty) * cols) + (r * regionWidth + j * T + t);
                            int dstIndex = ((Ty * (i * Ry + ry) + t) * cols) + (T * (j * R + r) + t);
                            output[dstIndex] = input[srcIndex];
                        }
                    }
                }
            }
        }
    } */
}


void openBinaryFile(const std::string& filename, cv::Mat& img, int rows, int cols) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo binario." << std::endl;
        exit(-1);
    }
    uint16_t* buffer = new uint16_t[rows * cols];
    file.read(reinterpret_cast<char*>(buffer), rows * cols * sizeof(uint16_t));
    if (file.gcount() != rows * cols * sizeof(uint16_t)) {
        std::cerr << "Error al leer los datos del archivo binario." << std::endl;
        delete[] buffer;
        exit(-1);
    }
    img = cv::Mat(rows, cols, CV_16UC1, buffer);
}

void deconstructAndReconstructImage(const cv::Mat& img, uchar* outputArray, int rows, int cols, const std::string& tapType) {
    uchar* inputArray = img.data;
    memset(outputArray, 0, rows * cols);
    applyTap(inputArray, outputArray, rows, cols, tapType);
}

int main() {

    cv::Mat img;
    int rows = 480, cols = 640;
    const std::string tapType = "4X2-1Y";   
    openBinaryFile("C:/CODE/VideoTaps/src/" + tapType + ".bin", img, rows, cols);

    double minVal, maxVal;
    cv::minMaxLoc(img, &minVal, &maxVal);
    cv::Mat img_normalizada;
    double alpha = 255.0 / (maxVal - minVal); // Contraste
    double beta = 50;  // Brillo (ajusta este valor)
    img.convertTo(img_normalizada, CV_8UC1, alpha, -minVal * alpha + beta);
    

    uchar* outputArray = new uchar[rows * cols];
    deconstructAndReconstructImage(img_normalizada, outputArray, rows, cols, tapType);

    cv::Mat reconstructed(rows, cols, CV_8UC1, outputArray);

    cv::imshow("Imagen Original", img_normalizada);
    cv::imshow("Imagen Reconstruida", reconstructed);
    cv::waitKey(0);

    cv::imwrite((tapType + "_Original.png"), img_normalizada);
    cv::imwrite((tapType + "_Reconstruida.png"), reconstructed);

    delete[] outputArray;


    //vector<string> testTaps = {"2X2", "4X", "2X2-1Y4", "2X2E-1Y2M"};

    auto [R, T, L, R2, T2, L2] = readTap(tapType);
    cout << "Tap: " << tapType << " -> R: " << R << ", T: " << T << ", L: " << L
         << " | R2: " << R2 << ", T2: " << T2 << ", L2: " << L2 << endl;


    return 0;
}
