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


void applyDDR(const uint8_t* input, uint16_t* out16, int rows, int cols, const string& tapType){
    auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);
    int simPixels = R * T * Ry * Ty;
    size_t cycles = size_t(rows) * cols / simPixels;
    
    std::vector<uint8_t> even(simPixels), odd(simPixels);
    const uint8_t* p = input;

    // Tenemos que permutar cada ciclo por el orden de llegada de los píxeles
    std::vector<int> perm(simPixels);
    for(int i = 0; i < simPixels; i++) perm[i] = i; // Permutación identidad

/*        
   // 1) 1X2-1Y2  → ciclo llega como [p10, p00, p11, p01] // TAP quiere [p00,p01,p10,p11]
    if (R==1 && T==2 && Ry==1 && Ty==2 && simPixels==4) { int m[4]={2,0,3,1}; for(int i=0;i<4;++i) perm[i]=m[i]; }
    // 2) 1X-2Y    → ciclo llega como [up, down], // TAP quiere [up, down]
    if (R==1 && T==1 && Ry==1 && Ty==2 && simPixels==2) { int m[2]={1,0};     for(int i=0;i<2;++i) perm[i]=m[i]; }
    // 3) 4XR-1Y   → ciclo llega como [1,0,3,2] (swap en cada pareja de R/2)
    if (R==4 && T==1 && Ry==1 && Ty==1 && simPixels==4){ int m[4]={1,0,3,2};  for(int i=0;i<4;++i) perm[i]=m[i]; }
*/
    for(int n = 0; n < cycles; n++){      
        // ---- CICLO 1: bits PARES ----
        std::memcpy(even.data(), p, simPixels); p += simPixels;
        // ---- CICLO 2: bits IMPARES ----
        std::memcpy(odd.data(), p, simPixels);  p += simPixels;

        
        for(int i = 0; i < simPixels; i++){
            int k = perm[i]; // Ajustamos orden de llegada
            out16[n * simPixels + i] = reconstruirPixel(even[k], odd[k]);
        }
    }

}

// VideoTap Standard (sin DDR)
template<typename T>
void applyTap(const T* input, T* output, int rows, int cols, const string& tapType) {
    auto [Rx, Tx, Lx, Ry, Ty, Ly] = readTap(tapType);

    const int regionWidth = cols / Rx;
    const int tapsX = regionWidth / Tx; 
    const int regionHeight = rows / Ry;
    const int tapsY = regionHeight / Ty;
    
    // buffers en T
    std::unique_ptr<T[]> bufferX(new T[rows * cols]);
    std::unique_ptr<T[]> bufferY(new T[rows * cols]);
    std::unique_ptr<T[]> buffer (new T[rows * cols]);

    std::memcpy(bufferX.get(), input,  size_t(rows)*cols*sizeof(T));
    std::memcpy(bufferY.get(), input,  size_t(rows)*cols*sizeof(T));
    std::memcpy(buffer .get(), input,  size_t(rows)*cols*sizeof(T));
    
    std::ofstream ofile1("step1.txt");
    std::ofstream ofile2("step2.txt");

    if(Rx != 1 || Tx != 1 || Lx != '\0'){
        if(Lx != '\0'){
            for(int i = 0; i < regionHeight; i++){
                for(int r = 0; r < Rx; r++){
                    for(int j = 0; j < tapsX; j++){
                        for(int t = 0; t < Tx; t++){
                            int srcIndex = (i * cols) + (r * regionWidth + j * Tx + t);
                            if(Lx == 'E'){
                                if(r >= Rx/2){
                                    if(Rx == 2){
                                        int dstIndex = (i * cols) + (regionWidth * (r + 1) - (Tx * j + t + 1));    
                                        bufferX[dstIndex] = input[srcIndex];
                                    } else {
                                        int dstIndex = (i * cols) + (regionWidth * (r + 1) - (Tx * (j + 1) - t));
                                        bufferX[dstIndex] = input[srcIndex];                                
                                    }
                                } else {
                                    int dstIndex = (i * cols) + r * regionWidth + j * Tx + t;
                                    bufferX[dstIndex] = input[srcIndex];
                                }
                            } else if (Lx == 'M'){
                                if(r >= Rx/2){
                                    int dstIndex = (i * cols) + r * regionWidth + j * Tx + t;
                                    bufferX[dstIndex] = input[srcIndex];
                                } else {
                                    int dstIndex = (i * cols) + (regionWidth * (r + 1) - (Tx * (j + 1) - t));
                                    bufferX[dstIndex] = input[srcIndex];
                                }
                            } else if (Lx == 'R'){
                                int dstIndex = (i * cols) + (regionWidth * (r + 1) - (Tx * (j + 1) - t));
                                bufferX[dstIndex] = input[srcIndex];

                            }
                        }
                    }
                }
            }
            // Lx = '\0';
        }
    
        for (int i = 0; i < regionHeight; i++){
            for (int r = 0; r < Rx; r++){
                for (int j = 0; j < tapsX; j++){
                    for (int t = 0; t < Tx; t++){
                        int srcIndex = (i * cols) + (r * regionWidth + j * Tx + t);
                        int dstIndex = (i * cols) + (Tx * (j * Rx + r) + t);
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
            std::memcpy(output, bufferY.get(), size_t(rows) * cols * sizeof(T));
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
            std::memcpy(output, bufferY.get(), size_t(rows) * cols * sizeof(T));
        }
    } else { 
        std::memcpy(output, buffer.get(), size_t(rows)*cols*sizeof(T));  
    } 
        
    ofile2.close();
    ofile1.close();
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
        std::vector<uint16_t> linear(N); // orden de llegada
        applyDDR(img.ptr<uint8_t>(), linear.data(), rows, cols, tapType);
        applyTap<uint16_t>(linear.data(), out16, rows, cols, tapType);
    } else {
        // Esperamos imagen 16b en orden "lineal" (o del tap de origen)
        CV_Assert(img.type() == CV_16UC1 && img.total() == N);
        applyTap<uint16_t>(img.ptr<uint16_t>(), out16, rows, cols, tapType);
    }
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

        const size_t N = size_t(rows) * cols;
        const size_t bytes = 2*N;   // (pares+impares)
        cv::Mat rawBytes(1, int(bytes), CV_16UC1);  
        std::ifstream f(path, std::ios::binary);
        f.read(reinterpret_cast<char*>(img8.data), std::streamsize(bytes));
        if (!f) { std::cerr << "Read failed\n"; return 1; }

        std::vector<uint16_t> out16(N);
        processFrame(rawBytes, out16.data(), rows, cols, tapType, isDDR);

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
        return 0;
        
        rows = 480, cols = 640; // .bin
        path = "C:/CODE/VideoTaps/src/input/" + tapType + ".bin";
        std::cout << "Rows: " << rows << "\nCols: " << cols << "\nPath: " << path << std::endl;
        // Abrimos
        cv::Mat img16;
        openBinaryFile("C:/CODE/VideoTap_Refactor/src/input/" + tapType + ".bin", img16, rows, cols);
        
        std::vector<uint16_t> out16(size_t(rows)*cols);
        processFrame(img16, out16.data(), rows, cols, tapType, isDDR);

        cv::Mat recon16(rows, cols, CV_16UC1, out16.data());

        saveAndShow(img8, recon16, tapType);
    }


    return 0;
}
