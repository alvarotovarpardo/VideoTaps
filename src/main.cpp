#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <opencv2/opencv.hpp>
#include <cstring>
#include <fstream>

struct TapInfo {
    char type; // 'X' o 'Y'
    int R;
    int T;
    char L;
};

std::vector<TapInfo> readTap(const std::string& tapType) {
    std::vector<TapInfo> taps;
    size_t pos = 0;
    while (pos < tapType.length()) {
        char type = tapType[pos + 1];
        if (type != 'X' && type != 'Y') {
            std::cerr << "Invalid format" << std::endl;
            return {};
        }
        int R = std::stoi(tapType.substr(pos, tapType.find(type, pos) - pos));
        pos = tapType.find(type, pos) + 1;
        int T = 1;
        if (pos < tapType.length() && std::isdigit(tapType[pos])) {
            size_t endT = pos;
            while (endT < tapType.length() && std::isdigit(tapType[endT])) {
                endT++;
            }
            T = std::stoi(tapType.substr(pos, endT - pos));
            pos = endT;
        }
        char L = '\0';
        if (pos < tapType.length() && std::isalpha(tapType[pos])) {
            L = tapType[pos];
            pos++;
        }
        if (pos < tapType.length() && tapType[pos] == '-') {
            pos++;
        }
        taps.push_back({type, R, T, L});
    }
    
    // Imprimir el resultado del parseo
    std::cout << "Parsing: " << tapType << std::endl;
    for (const auto& tap : taps) {
        std::cout << "R: " << tap.R << ", T: " << tap.T << ", L: " << (tap.L ? std::string(1, tap.L) : "None") << ", Axis: " << tap.type << std::endl;
    }
    std::cout << "-------------------" << std::endl;
    
    return taps;
}

void applyTap(const uchar* input, uchar* output, int rows, int cols, const std::string& tapType) {
    auto taps = readTap(tapType);
    std::memset(output, 0, rows * cols);
    
    for (const auto& tap : taps) {
        int regionSize = (tap.type == 'X') ? cols / tap.R : rows / tap.R;
        
        for (int r = 0; r < tap.R; r++) {
            for (int i = 0; i < regionSize; i++) {
                int srcIndex, dstIndex;
                if (tap.type == 'X') {
                    srcIndex = i * cols + r * regionSize;
                    dstIndex = i * cols + r;
                } else {
                    srcIndex = r * regionSize * cols + i;
                    dstIndex = i * cols + r * regionSize;
                }
                if (srcIndex < rows * cols && dstIndex < rows * cols) {
                    output[dstIndex] = input[srcIndex];
                }
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

int main() {
    cv::Mat img;
    int rows = 480, cols = 640;
    openBinaryFile("C:/CODE/VideoTaps/src/prueba.bin", img, rows, cols);

    double minVal, maxVal;
    cv::minMaxLoc(img, &minVal, &maxVal);
    cv::Mat img_normalizada;
    img.convertTo(img_normalizada, CV_8UC1, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));

    uchar* outputArray = new uchar[rows * cols];
    const std::string tapType = "2X2";
    applyTap(img_normalizada.data, outputArray, rows, cols, tapType);

    cv::Mat reconstructed(rows, cols, CV_8UC1, outputArray);

    cv::imshow("Imagen Original", img_normalizada);
    cv::imshow("Imagen Reconstruida", reconstructed);
    cv::waitKey(0);

    delete[] outputArray;
    return 0;
}
