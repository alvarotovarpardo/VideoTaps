#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <opencv2/opencv.hpp>
#include <cstring>
#include <fstream>
using namespace cv;


uint16_t reconstruirPixel(uint8_t pares, uint8_t impares) {
    uint16_t pixel = 0;
    for (int i = 0; i < 8; ++i) {
        pixel |= ((pares >> i) & 0x01) << (2 * i);
        pixel |= ((impares >> i) & 0x01) << (2 * i + 1);
    }
    return pixel;
}

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



void applyDDR(std::istream& input, std::vector<uint16_t>& output, int rows, int cols, const string& tapType){
	auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);
	int simPixels = R*T*Ry*Ty; // Pixeles simultáneos entrantes
    output.resize(rows*cols); // Ajustamos el output al tamaño del frame
	size_t cycles = rows*cols / simPixels; // Ciclos 'enviados'

	std::vector<uint8_t> even(simPixels), odd(simPixels);

	// For each cycle...
	for(int n = 0; n < cycles; n++){
		// CICLO 1: Píxeles pares
		input.read(reinterpret_cast<char*>(even.data()), simPixels); 
		// CICLO 2: Píxeles impares
		input.read(reinterpret_cast<char*>(odd.data()), simPixels);
		for(int i = 0; i < simPixels; i++){
			output[n * simPixels + i] = reconstruirPixel(even[i], odd[i]);
		}
	}
}


// VideoTap Standard (sin DDR)
void applyVideoTap(uchar* input, uchar* output, int rows, int cols, const string& tapType) {
    auto [R, T, L, Ry, Ty, Ly] = readTap(tapType);

    int regionWidth = cols / R;
    int tapsNumber = regionWidth / T; 
    int regionHeight = rows / Ry;
    int tapsNumberY = regionHeight / Ty;
    
    uchar* buffer = new uchar[rows * cols];
    memcpy(buffer, input, rows * cols);
    
    if(L != '\0'){
        // std::ofstream ofile("E2.txt");
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
                                   // ofile << srcIndex << " "<< dstIndex << "\n";// \t||\t " << r << " " << j << " " << t << " \n"; // << static_cast<int>(input[srcIndex]) << std::endl;
                                }
                            } else {
                                int dstIndex = (i * cols) + r * regionWidth + j * T + t;
                                buffer[dstIndex] = input[srcIndex];
                                // ofile << srcIndex << " "<< dstIndex << "\n"; // \t||\t " << r << " " << j <<  " " << t << " \n"; // << static_cast<int>(input[srcIndex]) << std::endl;
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
            } // ofile.close();
        }
    L = '\0';
    }
    



    if(Ty != 1 && T == 1){
        // std::ofstream ofile("Ly.txt");
        for (int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                int srcIndex = (i * cols) + j;
                int dstIndex = ((i - (i % 2)) * cols) + (2 * j + (i % 2));
                // ofile << srcIndex << " " << dstIndex << "\n";// \t||\t " << r << " " << j << " " << t << std::endl;
                output[dstIndex] = buffer[srcIndex];
            }
        } // ofile.close(); 
    } if(Ty != 1 && T != 1){
        // std::ofstream ofile("L_Ly.txt");
        for (int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                int srcIndex = (i * cols) + j;
                int dstIndex = ((i - (i % 2)) * cols) + (2 * j + (i % 2));
                // ofile << srcIndex << " " << dstIndex << "\n";// \t||\t " << r << " " << j << " " << t << std::endl;
                output[dstIndex] = buffer[srcIndex];
            }
        } // ofile.close(); 
    } else {
        // std::ofstream ofile("4X2E.txt");
        for (int i = 0; i < regionHeight; i++){
            for (int r = 0; r < R; r++){
                for (int j = 0; j < tapsNumber; j++){
                    for (int t = 0; t < T; t++){
                        int srcIndex = (i * cols) + (r * regionWidth + j * T + t);
                        int dstIndex = (i * cols) + (T * (j * R + r) + t);
                        output[dstIndex] = buffer[srcIndex];
                       // ofile << srcIndex << " " << dstIndex << "\n";// \t||\t " << r << " " << j << " " << t << std::endl;
                    }
                }
            } // ofile.close();
        }
    }
    delete [] buffer; buffer = nullptr;
}



int main()
{
    std:string tapType = "1X-1Y";
    bool hasDDR = false;
    int height;
    int width;




    return 0;
}