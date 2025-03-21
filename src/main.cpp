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
    
    cout << "[R, T, L] = [" << R << ", " << T << ", " << (L ? string(1, L) : "None") << "]\n";
    cout << "applying taps...\n";

    
    int regionWidth = cols / R;
    int numPackets = regionWidth / T; 
    
    for (int row = 0; row < rows; row++){
        uchar* newRow = new uchar[cols];
        memset(newRow, 0, cols);
        int idx = 0;
        
        
        for (int r = 0; r < R; r++){ // Regions
            int firstPixRegion = r * regionWidth;
    
            for (int p = 0; p < numPackets; p++){ // Packets
                for (int t = 0; t < T; t++){ // Pixels
                    int srcCol = firstPixRegion + p * T + t; 
                    newRow[idx++] = input[row * cols + srcCol]; 
                }
            }
        }
        memcpy(output + row * cols, newRow, cols * sizeof(uchar));
        delete [] newRow;
    }

}


int main() {
    // Abrir el archivo binario
    std::ifstream file("C:/CODE/VideoTaps/src/prueba.bin", std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "No se pudo abrir el archivo binario." << std::endl;
        return -1;
    }

    // Determina el tamaño de la imagen
    int rows = 480;  // Número de filas
    int cols = 640;  // Número de columnas

    // Crear un buffer para almacenar los datos del archivo binario
    uint16_t* buffer = new uint16_t[rows * cols]; // Buffer para valores de 16 bits

    // Leer el contenido del archivo binario
    file.read(reinterpret_cast<char*>(buffer), rows * cols * sizeof(uint16_t));
    if (file.gcount() != rows * cols * sizeof(uint16_t)) {
        std::cerr << "Error al leer los datos del archivo binario." << std::endl;
        delete[] buffer;
        return -1;
    }

    // Convertir el buffer en una matriz de OpenCV
    cv::Mat img(rows, cols, CV_16UC1, buffer); // CV_16UC1 es una imagen de 16 bits, un solo canal

    // Encontrar el valor mínimo y máximo en la imagen
    double minVal, maxVal;
    cv::minMaxLoc(img, &minVal, &maxVal);

    // Normalizar la imagen al rango de 0-255 (8 bits)
    cv::Mat img_normalizada;
    img.convertTo(img_normalizada, CV_8UC1, 255.0 / (maxVal - minVal), -minVal * 255.0 / (maxVal - minVal));

    // Mostrar la imagen ajustada
    cv::imshow("Imagen Ajustada", img_normalizada);
    cv::waitKey(0);

    // Liberar memoria
    delete[] buffer;
    return 0;
}
