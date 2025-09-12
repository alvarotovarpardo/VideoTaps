#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <opencv2/opencv.hpp>
#include <cstring>
#include <fstream>
using namespace cv;
using namespace std;

tuple<int, int, char, int, int, char> readTap(const string& tapType) { // 4XR-1Y // 2X-1Y // 2X-1YE // 1X
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

    return {R, T, L, R2, T2, L2}; // RXTL - R2YT2L2
                                  // 4X2E - 2Y2M
}

uint16_t reconstruirPixel(uint8_t pares, uint8_t impares) {
    uint16_t pixel = 0;
    for (int i = 0; i < 8; ++i) {
        pixel |= ((pares >> i) & 0x01) << (2 * i);
        pixel |= ((impares >> i) & 0x01) << (2 * i + 1);
    }
    return pixel;
}

void openBinaryFile(const std::string& filename, cv::Mat& img, int rows, int cols) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("Error opening " + filename);
    const size_t need = size_t(rows)*cols*sizeof(uint16_t);
    if (size_t(file.tellg()) < need) throw std::runtime_error("Archivo mas pequeño que un frame");
    file.seekg(0);

    img.create(rows, cols, CV_16UC1);
    file.read(reinterpret_cast<char*>(img.data), need);
    if (size_t(file.gcount()) != need) throw std::runtime_error("Lectura incompleta");
}


void applyDDR(uchar* input, uchar* output, int rows, int cols, const string& tapType){
    auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);
    int simPixels = R*T*Ry*Ty;      // Pixeles simultáneos por ciclo
    //output.resize(rows*cols);     // Ajustamos el output al tamaño del frame
    size_t cycles = rows*cols / simPixels; // Ciclos 'enviados'

    uint8_t* p = reinterpret_cast<uint8_t*>(input);

    // For each cycle...
    for(int n = 0; n < cycles; n++){
        // CICLO 1: Píxeles pares
        //input.read(reinterpret_cast<char*>(even.data()), simPixels); 
        uint8_t* even = p;
        // CICLO 2: Píxeles impares
        uint8_t* odd = p + simPixels;
        //input.read(reinterpret_cast<char*>(odd.data()), simPixels);
        for(int i = 0; i < simPixels; i++){
            uint16_t pix = reconstruirPixel(even[i], odd[i]);
            output[n * simPixels + i] = static_cast<uchar>(pix); // (p >> 8)
        }
        p+= 2 * simPixels; // Siguiente ciclo (factor 2 : par/impar)
    }
}

// VideoTap Standard (sin DDR)
void applyTap(uchar* input, uchar* output, int rows, int cols, const string& tapType) {
    auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);

    int regionWidth = cols / R;
    int tapsNumber = regionWidth / T; 
    int regionHeight = rows / Ry;
    int tapsNumberY = regionHeight / Ty;
    
    uchar* buffer = new uchar[rows * cols];
    memcpy(buffer, input, rows * cols);
    
    if(L != '\0'){
        for(int i = 0; i < regionHeight; i++){
            for(int r = 0; r < R; r++){
                for(int j = 0; j < tapsNumber; j++){
                    for(int t = 0; t < T; t++){
                        int srcIndex = (i * cols) + (r * regionWidth + j * T + t);
                        if(L == 'E'){
                            if(r >= R/2){
                                if(R == 2){
                                    int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * j + t + 1));    
                                    buffer[dstIndex] = input[srcIndex];
                                } else {
                                    int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * (j + 1) - t));
                                    buffer[dstIndex] = input[srcIndex];                                
                                }
                            } else {
                                int dstIndex = (i * cols) + r * regionWidth + j * T + t;
                                buffer[dstIndex] = input[srcIndex];
                            }
                        } else if (L == 'M'){
                            if(r >= R/2){
                                int dstIndex = (i * cols) + r * regionWidth + j * T + t;
                                buffer[dstIndex] = input[srcIndex];
                            } else {
                                int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * (j + 1) - t));
                                buffer[dstIndex] = input[srcIndex];
                            }
                        } else if (L == 'R'){
                            int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * (j + 1) - t));
                            buffer[dstIndex] = input[srcIndex];

                        }
                    }
                }
            }
        }
    L = '\0';
    }
    



    if(Ty != 1 && T == 1){
        for (int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                int srcIndex = (i * cols) + j;
                int dstIndex = ((i - (i % 2)) * cols) + (2 * j + (i % 2));
                output[dstIndex] = buffer[srcIndex];
            }
        } 
    } else if(Ty != 1 && T != 1){
        for (int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                int srcIndex = (i * cols) + j;
                int dstIndex = ((i - (i % 2)) * cols) + (2 * j + (i % 2));
                output[dstIndex] = buffer[srcIndex];
            }
        }
    } else {
        for (int i = 0; i < regionHeight; i++){
            for (int r = 0; r < R; r++){
                for (int j = 0; j < tapsNumber; j++){
                    for (int t = 0; t < T; t++){
                        int srcIndex = (i * cols) + (r * regionWidth + j * T + t);
                        int dstIndex = (i * cols) + (T * (j * R + r) + t);
                        output[dstIndex] = buffer[srcIndex];
                    }
                }
            } 
        }
    } 
    delete [] buffer; buffer = nullptr;
}

cv::Mat normalizeImage(const cv::Mat& img) {
    double minVal, maxVal; cv::minMaxLoc(img, &minVal, &maxVal);
    cv::Mat out;

    if (maxVal > minVal) {
        double alpha = 65535.0 / (maxVal - minVal);
        double beta  = -minVal * alpha;          // sin offset extra
        img.convertTo(out, CV_16UC1, alpha, beta);
    } else {
        out.setTo(img);                             // imagen constante
    }
    return out;
}


void processFrame(const cv::Mat& img, uchar* outputArray, int rows, int cols, const std::string& tapType, bool isDDR) {
    uchar* inputArray = img.data;
    const size_t N = size_t(rows) * cols;
    // buffer intermedio
    std::unique_ptr<uchar[]> tmp(new uchar[N]);
    memset(tmp.get(), 0, N);

    if(isDDR){
        applyDDR(inputArray, tmp.get(), rows, cols, tapType);
    } else {
        std::memcpy(tmp.get(), inputArray, N);
    }

    applyTap(tmp.get(), outputArray, rows, cols, tapType);
}

void saveAndShow(cv::Mat& input, cv::Mat& output, const std::string& tapType){
    cv::imshow("Imagen Original", input);
    cv::imshow("Imagen Reconstruida", output);
    cv::waitKey(0);

    //cv::imwrite((tapType + "_Original.png"), input);
    //cv::imwrite((tapType + "_Reconstruida.png"), output);
}

int main() {

    cv::Mat img, img_norm16;
    int rows = 480, cols = 640; // .bin
    
    size_t frameSize = rows * cols * 2;


    std::string tapType; 
    std::cout << "Introduce Tap Geometry: "; std::getline(std::cin, tapType);
    
    openBinaryFile("C:/CODE/VideoTap_Refactor/src/input/" + tapType + ".bin", img, rows, cols);

    img_norm16 = normalizeImage(img); // CV_16UC1 -> CV_8UC1
    cv::Mat img8; img_norm16.convertTo(img8, CV_8UC1, 1.0/256.0);

    std::unique_ptr<uchar[]> out(new uchar[rows*cols]);
    processFrame(img8, out.get(), rows, cols, tapType, /*isDDR=*/false);

    cv::Mat reconstructed(rows, cols, CV_8UC1, out.get());
    saveAndShow(img8, reconstructed, tapType);

/*
    // Guardar la primera fila de img_normalizada en un archivo .txt
    std::ofstream outFile(tapType + "_primera_fila.txt");
    if (outFile.is_open()) {
        for (int col = 0; col < cols; ++col) {
            outFile << static_cast<int>(img_normalizada.at<uchar>(0, col));
            if (col < cols - 1) outFile << " ";
        }
        outFile.close();
        std::cout << "Primera fila guardada en " << tapType + "_primera_fila.txt" << std::endl;
    } else {
        std::cerr << "No se pudo abrir el archivo para escritura." << std::endl;
    }


    //vector<string> testTaps = {"2X2", "4X", "2X2-1Y4", "2X2E-1Y2M"};

    auto [R, T, L, R2, T2, L2] = readTap(tapType);
    cout << "Tap: " << tapType << " -> R: " << R << ", T: " << T << ", L: " << L
         << " | R2: " << R2 << ", T2: " << T2 << ", L2: " << L2 << endl;
*/
    
    return 0;
}
