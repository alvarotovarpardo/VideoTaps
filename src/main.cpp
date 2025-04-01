#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstring> 
#include <string>
#include <sstream>
#include <fstream>
#include <tuple>
using namespace cv;
using namespace std;

tuple<int, int, char> readTap(const string& tapType) {
    size_t xPos = tapType.find('X');
    if(xPos == string::npos || xPos == 0){
        cerr << "Invalid format" << endl;
        return {-1,-1,'\0'};
    }

    int R = stoi(tapType.substr(0,xPos));
    int T = 1;
    char L = '\0';

    if(xPos + 1 < tapType.length()){
        size_t posL = xPos + 1;
        while (posL < tapType.length() && isdigit(tapType[posL])) {
            posL++;
        }

        if (posL > xPos + 1) {
            T = stoi(tapType.substr(xPos + 1, posL - (xPos + 1)));
        }

        if (posL < tapType.length()) {
            L = tapType[posL];
        }
    }
    
    return {R, T, L};
}

void applyTap(const uchar* input, uchar* output, int rows, int cols, const string& tapType){
    auto [R, T, L] = readTap(tapType);
    
    int regionWidth = cols / R; 
    
    for (int row = 0; row < rows; row++) {
        for (int r = 0; r < R; r++) { // Región por Tap
            for (int col = 0; col < regionWidth; col++) { 
                int srcIndex = row * cols + r * regionWidth + col;
                int dstIndex = row * cols + col * R + r; 
                output[dstIndex] = input[srcIndex];
            }
        }
    }
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
    openBinaryFile("C:/CODE/VideoTaps/src/4X.bin", img, rows, cols);

    double minVal, maxVal;
    cv::minMaxLoc(img, &minVal, &maxVal);
    cv::Mat img_normalizada;
    img.convertTo(img_normalizada, CV_8UC1, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));

    uchar* outputArray = new uchar[rows * cols];
    const std::string tapType = "4X";
    deconstructAndReconstructImage(img_normalizada, outputArray, rows, cols, tapType);

    cv::Mat reconstructed(rows, cols, CV_8UC1, outputArray);

    cv::imshow("Imagen Original", img_normalizada);
    cv::imshow("Imagen Reconstruida", reconstructed);
    cv::waitKey(0);

    delete[] outputArray;
    return 0;
}