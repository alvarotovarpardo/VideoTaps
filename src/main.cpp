#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstring> 
#include <string>
#include <sstream>
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
    Mat img = imread("C:/CODE/VideoTaps/src/image.jpg", IMREAD_GRAYSCALE);
    if (img.empty()) {
        std::cerr << "Error al cargar la imagen." << std::endl;
        return -1;
    }

    int rows = img.rows;
    int cols = img.cols;

    uchar* inputArray = img.data;
    uchar* outputArray = new uchar[rows * cols];
    memset(outputArray, 0, rows * cols);
    
    const string tapType = "2X2"; 
    
    applyTap(inputArray, outputArray, rows, cols, tapType);
    
    Mat reconstructed(rows, cols, CV_8UC1, outputArray);

    imshow("Imagen Original", img);
    imshow("Imagen Reconstruida", reconstructed);
    
    waitKey(0);
    return 0;
}