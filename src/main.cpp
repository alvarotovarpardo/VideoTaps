#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
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

void DDR_1X_2Y(const uint8_t* input, uint16_t* out16, int rows, int cols){
        
    std::vector<uint8_t> even(2), odd(2);
    const uint8_t* p = input;

    for (int i = 0; i < rows; i += 2) {
        for (int j = 0; j < cols; j++) {
            std::memcpy(even.data(), p, 2); p += 2;
            std::memcpy(odd.data(),  p, 2); p += 2;

            uint16_t pUp = reconstruirPixel(even[0], odd[0]); 
            uint16_t pDown = reconstruirPixel(even[1], odd[1]); 

            out16[i*cols + j]         = pDown;
            out16[(i+1)*cols + j]     = pUp;
        }
    }
}

void DDR_1X2_1Y2(const uint8_t* input, uint16_t* out16, int rows, int cols){
        
    std::vector<uint8_t> even(4), odd(4);
    const uint8_t* p = input;

    for (int i = 0; i < rows; i += 2) {
        for (int j = 0; j < cols; j += 2) {
            std::memcpy(even.data(), p, 4); p += 4;
            std::memcpy(odd.data(),  p, 4); p += 4;
            
            // (im)pares[0] → (1,0), (im)pares[1] → (0,0), (im)pares[2] → (1,1), (im)pares[3] → (0,1)
            uint16_t p10 = reconstruirPixel(even[0], odd[0]);
            uint16_t p00 = reconstruirPixel(even[1], odd[1]);
            uint16_t p11 = reconstruirPixel(even[2], odd[2]);
            uint16_t p01 = reconstruirPixel(even[3], odd[3]);

            out16[i      * cols + (j    )] = p00;  // (0,0)
            out16[(i + 1)* cols + (j    )] = p10;  // (1,0)
            out16[i      * cols + (j + 1)] = p01;  // (0,1)
            out16[(i + 1)* cols + (j + 1)] = p11;  // (1,1)
        }
    }
}

void DDR_4XR_1Y(const uint8_t* input, uint16_t* out16, int rows, int cols){
        
    std::vector<uint8_t> even(4), odd(4);
    const uint8_t* p = input;
    
    int regionWidth = cols / 4;       // 80   
    const int regionOrder[4] = {2, 1, 4, 3}; // Orden de llegada: 159, 79, 319, 239, 158, 78, 318, 238, ...

    for (int i = 0; i < rows; i++) {
        // cada j agrupa 4 píxeles (uno por cada canal/región)
        for (int j = 0; j < regionWidth; j++) {
            std::memcpy(even.data(), p, 4); p += 4;
            std::memcpy(odd.data(),  p, 4); p += 4;
            
            for (int r = 0; r < 4; ++r) {
                int region = regionOrder[r];
                int col = region * regionWidth - j - 1;
                size_t idx = i * cols + col;    // posición en el array lineal
                out16[idx] = reconstruirPixel(even[r], odd[r]);
            }
        }
    }
}

void applyDDR(const uint8_t* input, uint16_t* out16, int rows, int cols, const std::string& tapType)
{
    auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);
    const int simPixels = R * T * Ry * Ty;          // 4 en 1X2-1Y2
    const size_t N = size_t(rows) * cols;
    const size_t cycles = N / simPixels;


    if (R==1 && T==2 && Ry==1 && Ty==2) {DDR_1X2_1Y2(input, out16, rows, cols); return;}
    if (R==1 && T==1 && Ry==2 && Ty==1) {DDR_1X_2Y(input, out16, rows, cols); return;}
    if (R==4 && T==1 && Ry==1 && Ty==1) {DDR_4XR_1Y(input, out16, rows, cols); return;}


/*
    std::vector<uint8_t> even(simPixels), odd(simPixels);
    const uint8_t* p = input;

    // Para 1X2-1Y2 dejamos identidad: [p10,p00,p11,p01] tal cual.
    for (size_t n = 0; n < cycles; ++n) {
        std::memcpy(even.data(), p, simPixels); p += simPixels;
        std::memcpy(odd.data(),  p, simPixels); p += simPixels;
        for (int i = 0; i < simPixels; ++i) {
            out16[n * simPixels + i] = reconstruirPixel(even[i], odd[i]);
        }
    }
*/
}


// VideoTap Standard (sin DDR)
template<typename T>
void invertXtap(const T* input, T* output, int rows, int cols, int Rx, int Tx, int Lx){
    const int regionWidth = cols / Rx;
    const int tapsX = regionWidth / Tx;

    if (Lx == '\0') {
        std::memcpy(output, input, size_t(rows)*cols*sizeof(T));
        return;
    }

    auto invX = [&](int r)->bool{
        if (Lx == '\0') return false;
        const bool right = (r >= Rx/2), left = !right;
        if (Lx == 'R') return true;
        if (Lx == 'E') return right;
        if (Lx == 'M') return left;
        return false;
    };

    // Gestiona el envío reverso de bloques dentro de un tap Tx
    auto innerReverse = [&](int r)->bool{
        if (!(Lx == 'E' || Lx == 'R')) return false;
        bool needReverse = true;
        if (Lx == 'E' && Rx == 4 && Tx == 2) needReverse = false;
        return needReverse && invX(r);
    };

    for(int i = 0; i < rows; i++){
        const T* srcRow = input + size_t(i) * cols;
              T* dstRow = output + size_t(i) * cols;

        for(int rx = 0; rx < Rx; rx++){
            const bool invert = invX(rx); // Invert region
            const bool innvert = innerReverse(rx); // Invert pixel order

            const T* srcRegion = srcRow + rx * regionWidth;
                  T* dstRegion = dstRow + rx * regionWidth;

            if(!invert){
                std::memcpy(dstRegion, srcRegion, size_t(regionWidth)*sizeof(T));
                continue;
            }

            if(innvert){
                // Reinvertimos intra-bloque
                for(int j = 0; j < regionWidth; j++){
                    dstRegion[j] = srcRegion[regionWidth - j - 1];
                }
            } else {
                // Reinvertimos el orden de bloques de tamaño Tx
                for(int j = 0; j < tapsX; j++){
                    const int srcJ = (tapsX - 1 - j) * Tx;
                    std::memcpy(dstRegion + j*Tx, srcRegion + srcJ, size_t(Tx) * sizeof(T));
                }
            }
        }
    }
}

template<typename T>
void applyXtap(const T* input, T* output, int rows, int cols, int Rx, int Tx, char Lx){

    const int regionWidth = cols / Rx;
    const int tapsX = regionWidth / Tx;

    std::unique_ptr<T[]> tmp(new T[size_t(rows) * cols]);
    T* buffer = tmp.get();
    
    if(Lx != '\0'){
        invertXtap(input, buffer, rows, cols, Rx, Tx, Lx);
    } else {
        std::memcpy(buffer, input, size_t(rows)*cols*sizeof(T));
    }

    for(int i = 0; i < rows; i++){
        const T* srcRow = buffer + size_t(i) * cols;
              T* dstRow = output + size_t(i) * cols;

        for(int r = 0; r < Rx; r++){            
            const T* srcRegion = srcRow + r * regionWidth;

            for (int j = 0; j < tapsX; j++){
                const T* srcBlock = srcRegion + j * Tx;
                      T* dstBlock = dstRow + Tx * (j*Rx + r);
                
                for(int tx = 0; tx < Tx; tx++){
                    dstBlock[tx] = srcBlock[tx];
                }
            }
        }
    }  
}

template<typename T>
void invertYtap(const T* input, T* output, int rows, int cols, int Ry, char Ly) {
    const int regionHeight = rows / Ry;

    if (Ly == '\0') {
        std::memcpy(output, input, size_t(rows)*cols*sizeof(T));
        return;
    }

    auto invY = [&](int ry)->bool{
        if (Ly == '\0') return false;
        const bool bottom = (ry >= Ry/2), top = !bottom;
        if (Ly=='R') return true;
        if (Ly=='E') return bottom;
        if (Ly=='M') return top;
        return false;
    };
    for (int ry = 0; ry < Ry; ++ry) {
        const bool flip = invY(ry);
        for (int i = 0; i < regionHeight; ++i) {
            const int srcRow = flip
                ? (ry*regionHeight + (regionHeight - 1 - i))  // espejo vertical dentro de la región
                : (ry*regionHeight + i);
            const int dstRow = ry*regionHeight + i;
            std::memcpy(output + size_t(dstRow)*cols, input  + size_t(srcRow)*cols, size_t(cols)*sizeof(T));
        
        }
    }
}

template<typename T>
void applyYtap(const T* input, T* output, int rows, int cols, int Ry, int Ty, char Ly) {
    const int regionHeight = rows / Ry;
    
    std::unique_ptr<T[]> tmp(new T[size_t(rows) * cols]);
    T* buffer = tmp.get();
    
    if(Ly != '\0'){
        invertYtap(input, buffer, rows, cols, Ry, Ly);
    } else {
        std::memcpy(buffer, input, size_t(rows)*cols*sizeof(T));
    }

    for (int i = 0; i < regionHeight; ++i) {
        for (int ry = 0; ry < Ry; ++ry) {
            const int lane   = i % Ty; 
            const int srcRow = ry*regionHeight + i;
            const int dstRow = i * Ry - i % Ty; 

            const T* src = buffer  + size_t(srcRow)*cols;
                  T* dst = output + size_t(dstRow)*cols;

            for (int j = 0; j < cols; j++) {
                const int srcCol = j;
                const int dstCol = Ty*Ry*j + ry + lane;
                dst[dstCol] = src[srcCol];
            }
        }
    } 
}

template<typename T>
void applyMixTap(const T* input, T* output, int rows, int cols, int Ry, int Ty, int Rx, int Tx, char Ly) {
    const int regionHeight = rows / Ry;
    
    std::unique_ptr<T[]> tmp(new T[size_t(rows) * cols]);
    T* buffer = tmp.get();
    
    if(Ly != '\0'){
        invertYtap(input, buffer, rows, cols, Ry, Ly);
    } else {
        std::memcpy(buffer, input, size_t(rows)*cols*sizeof(T));
    }

    for (int i = 0; i < regionHeight; ++i) {
        for (int ry = 0; ry < Ry; ++ry) {
            // const int lane   = i % Ty; 
            const int srcRow = ry*regionHeight + i;
            const int dstRow = i * Ry - i % Ty; 

            const T* src = buffer  + size_t(srcRow)*cols;
                  T* dst = output + size_t(dstRow)*cols;

            for(int rx = 0; rx < Rx; rx++){
                for (int j = 0; j < cols/(Rx*Tx); j++) {
                    for(int tx = 0; tx < Tx; tx++){
                        const int srcCol = rx*cols/Rx + j*Tx + tx;
                        const int dstCol = (Ry*Rx*Tx*Ty)*j + (Rx*ry + rx)*Tx + i%Ty + tx;
                        dst[dstCol] = src[srcCol];
                        /*if(j < 9 && i < 2) std::cout << 
                            "(" << i << "," << ry << "," << rx << "," << j << "): " << 
                            srcCol << " -> " << dstCol << std::endl;
                        */
                    }
                }
            }
        }
    }
}

template<typename T>
void applyTap(const T* input, T* output, int rows, int cols, const string& tapType) {
    auto [Rx, Tx, Lx, Ry, Ty, Ly] = readTap(tapType);
    
    bool hasYtap = (Ry != 1 || Ty != 1 || Ly != '\0');
    bool hasXtap = (Rx != 1 || Tx != 1 || Lx != '\0');
    bool hasRyTx = (Ry != 1 && Tx != 1);
    bool hasXYtap = hasYtap && hasXtap;
    bool hasComplexTap = (Lx != '\0' || Ly != '\0') && hasXYtap;

    if(hasComplexTap){
        applyMixTap(input, output, rows, cols, Ry, Ty, Rx, Tx, Ly);
        return;
    }

    if(hasXYtap){            
        std::unique_ptr<T[]> tmp(new T[size_t(rows) * cols]);
        T* buffer = tmp.get();
        applyYtap(input, buffer, rows, cols, Ry, Ty, Ly);
        applyXtap(buffer, output, rows, cols, Rx, Tx, Lx);
        return;
    } else {
        if(hasYtap){
            applyYtap(input, output, rows, cols, Ry, Ty, Ly);
            return;
        }    
        applyXtap(input, output, rows, cols, Rx, Tx, Lx);
        return;
    }
    
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


void processFrame(const cv::Mat& img, uint16_t* out16, int rows, int cols, const std::string& tapType, bool isDDR)
{
    const size_t N = size_t(rows) * cols;

    if (isDDR) {
        // Esperamos 2*N bytes (pares+impares)
        CV_Assert(img.type() == CV_8UC1 && img.total() == 2*N);
        applyDDR(img.ptr<uint8_t>(), out16, rows, cols, tapType);  
        //applyTap<uint16_t>(linear.data(), out16, rows, cols, tapType);
    } else {
        CV_Assert(img.type() == CV_16UC1 && img.total() == N);
        applyTap<uint16_t>(img.ptr<uint16_t>(), out16, rows, cols, tapType);
    }
}

void saveAndShow(cv::Mat& input, cv::Mat& output, const std::string& tapType){
    cv::imshow("Imagen Original", input);
    cv::imshow("Imagen Reconstruida", output);
    cv::waitKey(0);
    // cv::imwrite(("Imagen_Original.png"), input);
    // cv::imwrite((tapType + ".png"), output);

}

int main() {

    cv::Mat img, img8, img_norm16;
    // DDR Videotap
    
    /*std::string askDDR; 
    std::cout << "Is DDR? (0/1): "; std::getline(std::cin, askDDR);
    bool isDDR = (askDDR == "1");
    */
    bool isDDR = false;

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

        const size_t N = size_t(rows)*cols;
        const size_t bytes = 2*N;
        cv::Mat rawBytes(1, int(bytes), CV_8UC1);
        std::ifstream f(path, std::ios::binary);
        f.read(reinterpret_cast<char*>(rawBytes.data), std::streamsize(bytes));
        if (!f) { std::cerr << "Read failed\n"; return 1; }

        std::vector<uint16_t> out16(N);
        
        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();

        processFrame(rawBytes, out16.data(), rows, cols, tapType, isDDR);
        auto t1 = clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        std::cout << "processFrame DDR time: " << ms << " ms\n";
        
        cv::Mat recon16(rows, cols, CV_16UC1, out16.data());

        // Normaliza
        double minV, maxV;
        cv::minMaxLoc(recon16, &minV, &maxV);
        cv::Mat recon8;
        if (maxV > minV) {
            recon16.convertTo(recon8, CV_8UC1, 255.0 / (maxV - minV), -minV * (255.0 / (maxV - minV)));
        } else {
            recon16.convertTo(recon8, CV_8UC1, 1.0/256.0);
        }

        cv::imshow("Imagen Reconstruida (8-bit)", recon8);
        cv::waitKey(0);            
        
    } else {        
        rows = 480, cols = 640; // .bin
        // Abrimos
        cv::Mat img16;
        openBinaryFile("C:/CODE/VideoTaps/src/input/" + tapType + ".bin", img16, rows, cols);
        
        std::vector<uint16_t> out16(size_t(rows)*cols);

        using clock = std::chrono::steady_clock;
        auto t0 = clock::now();
        
        processFrame(img16, out16.data(), rows, cols, tapType, isDDR);
        auto t1 = clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        std::cout << "processFrame time: " << ms << " ms\n";

        cv::Mat recon16(rows, cols, CV_16UC1, out16.data());
        cv::Mat recon16norm = normalizeImage(recon16); // Output normalizado
        cv::Mat img16norm   = normalizeImage(img16);   // Input normalizado
        saveAndShow(img16norm, recon16norm, tapType);
    }


    return 0;
}
