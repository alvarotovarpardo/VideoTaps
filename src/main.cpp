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
            output[n * simPixels + i] = static_cast<uchar>(pix >> 8); // (p >> 8)
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
    
    uchar* bufferX = new uchar[rows * cols];
    memcpy(bufferX, input, rows * cols);
    
    uchar* bufferY = new uchar[rows * cols];
    memcpy(bufferY, input, rows * cols);

    uchar* buffer = new uchar[rows * cols];
    memcpy(buffer, input, rows * cols);
    
    std::ofstream ofile1("step1.txt");
    std::ofstream ofile2("step2.txt");

    if(R != 1 || T != 1 || L != '\0'){
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
                                        bufferX[dstIndex] = input[srcIndex];
                                    } else {
                                        int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * (j + 1) - t));
                                        bufferX[dstIndex] = input[srcIndex];                                
                                    }
                                } else {
                                    int dstIndex = (i * cols) + r * regionWidth + j * T + t;
                                    bufferX[dstIndex] = input[srcIndex];
                                }
                            } else if (L == 'M'){
                                if(r >= R/2){
                                    int dstIndex = (i * cols) + r * regionWidth + j * T + t;
                                    bufferX[dstIndex] = input[srcIndex];
                                } else {
                                    int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * (j + 1) - t));
                                    bufferX[dstIndex] = input[srcIndex];
                                }
                            } else if (L == 'R'){
                                int dstIndex = (i * cols) + (regionWidth * (r + 1) - (T * (j + 1) - t));
                                bufferX[dstIndex] = input[srcIndex];

                            }
                        }
                    }
                }
            }
            L = '\0';
        }
    
        for (int i = 0; i < regionHeight; i++){
            for (int r = 0; r < R; r++){
                for (int j = 0; j < tapsNumber; j++){
                    for (int t = 0; t < T; t++){
                        int srcIndex = (i * cols) + (r * regionWidth + j * T + t);
                        int dstIndex = (i * cols) + (T * (j * R + r) + t);
                        buffer[dstIndex] = bufferX[srcIndex];
                        ofile1 << srcIndex << " " << dstIndex << "\n";
                    }
                }
            } 
        }
        ofile1.close();
    }

    if(Ry != 1 || Ty != 1){
        if(Ry != 1){
            size_t srcIndex = 0;
            for(int i = 0; i < rows; i+=2){ // Ry = 2 => mitad de filas
                for(int j = 0; j < cols; j++){
                    for(int ry = 0; ry < Ry; ry++){
                        int dstIndex = (i + ry) * cols + j; // TODO: parecido a 2X-1Y2, unificar
                        bufferY[srcIndex++] = input[dstIndex];
                    }
                }
            }
            memcpy(output, bufferY, rows * cols);
        }

        if(Ty != 1){ 
            for (int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    int srcIndex = (i * cols) + j;
                    int dstIndex = ((i - (i % 2)) * cols) + (2 * j + (i % 2));
                    bufferY[dstIndex] = buffer[srcIndex];
                    ofile2 << srcIndex << " " << dstIndex << "\n";
                }
            }
            memcpy(output, bufferY, rows * cols);
        }
    } else { 
        memcpy(output, buffer, rows * cols);    
    } 
        
    ofile2.close();
    
/*
    // 2X-1Y2
    size_t src = 0;
    if(Ty != 1){ 
    for (int i = 0; i < rows; i += Ty) {          
        for (int j = 0; j < regionWidth; ++j) {   
            for (int ty = 0; ty < Ty; ++ty) {     
                for (int r = 0; r < R; ++r) { 
                    int dstRow = i + ty;
                    int dstCol = r * regionWidth + j;
                    int dstIndex = (i + ty) * cols + r * regionWidth + j;
                    ofile1 << src << " " << dstIndex << "\n";
                    output[src++] = input[dstIndex];
                    
                    }
                } 
            }
        }
        memcpy(output, buffer, rows * cols);
    } else {
        memcpy(output, buffer, rows * cols);    
    } 
*/
    ofile1.close();

    delete [] buffer; buffer = nullptr;
    delete [] bufferX; bufferX = nullptr;
    delete [] bufferY; bufferY = nullptr;
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

    cv::Mat img, img8, img_norm16;
    // DDR Videotap
    
    std::string askDDR; 
    std::cout << "Is DDR? (0/1): "; std::getline(std::cin, askDDR);
    bool isDDR = (askDDR == "1");

    std::string tapType; 
    std::cout << "Introduce Tap Geometry: "; std::getline(std::cin, tapType);
    int rows, cols;
    std::string path;

    // Esto tiene que ver con la lectura de .raw, no debería afectar a la llegada de un frame en streaming
    if(isDDR){
        if(tapType == "4XR-1Y"){
            rows = 240, cols = 320; // .raw Global Shutter
        } else {
            rows = 1024, cols = 1280; // .raw
        }
        path = "C:/CODE/VideoTaps/src/input/0_DDR_" + tapType + ".raw";
        std::cout << "Rows: " << rows << "\nCols: " << cols << "\nPath: " << path << std::endl;
        size_t frameSize = rows * cols * 2;

        const size_t bytes = size_t(rows) * cols * 2;   // 2*N (pares+impares)
        img8.create(1, int(bytes), CV_8UC1);  

        std::ifstream f(path, std::ios::binary);
        f.read(reinterpret_cast<char*>(img8.data), std::streamsize(bytes));
        if (!f) { std::cerr << "Read failed\n"; return 1; }
        
    } else {
        rows = 480, cols = 640; // .bin
        path = "C:/CODE/VideoTaps/src/input/" + tapType + ".bin";
        std::cout << "Rows: " << rows << "\nCols: " << cols << "\nPath: " << path << std::endl;
        // Abrimos
        openBinaryFile("C:/CODE/VideoTap_Refactor/src/input/" + tapType + ".bin", img, rows, cols);
        // Normalizamos
        img_norm16 = normalizeImage(img); // CV_16UC1 -> CV_8UC1
        img_norm16.convertTo(img8, CV_8UC1, 1.0/256.0);

    }

    std::unique_ptr<uchar[]> out(new uchar[rows*cols]);
    processFrame(img8, out.get(), rows, cols, tapType, isDDR);

    cv::Mat reconstructed(rows, cols, CV_8UC1, out.get());
    saveAndShow(img8, reconstructed, tapType);
    
    return 0;
}
