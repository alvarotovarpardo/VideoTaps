#include <opencv2/opencv.hpp>
#include <iostream>
#include <cstring> 
#include <string>
#include <sstream>
using namespace cv;
using namespace std;

Mat readTap(const Mat& img, const string tapType) {
    
    // Searching for X to get R
    size_t xPos = tapType.find('X');
    if(xPos == std::string::npos || xPos == 0){
        std::cerr << "Invalid format" << std::endl;
        return cv::Mat();
    }

    // <R>X<T><L>

    // Default
    int R = 0;
    int T = 0;
    char L = '\0';

    R = std::stoi(tapType.substr(0,xPos));

    // Searching if T and L exists
    if(xPos + 1 < tapType.length()){
        size_t posL = xPos + 1;
        while (posL < tapType.length() && isdigit(tapType[posL])) {
            posL++;
        }

        // Extract T if exists
        if (posL > xPos + 1) {
            T = std::stoi(tapType.substr(xPos + 1, posL - (xPos + 1)));
        }

        // Extract L if exists
        if (posL < tapType.length()) {
            L = tapType[posL];
        }
    }
    


    std::cout << "[R, T, L] = [" << R << ", " << T << ", " << L << "]\n"; 



    Mat reconstructed = Mat::zeros(img.size(), img.type());
    
    cout << "Reconstruccion de la imagen completada." << endl;
    return reconstructed;
}

int main() {
    
    Mat img = imread("C:/CODE/SensiaNuc_STL/Test_Environment/src/image.jpg", IMREAD_GRAYSCALE);
    if (img.empty()) {
        cerr << "Error al cargar la imagen." << endl;
        return -1;
    }
    
    const std::string tapType = "10X40E";
    Mat reconstructed = readTap(img, tapType);
    
    // Mostrar la imagen original y la reconstruida para verificar que se procesó correctamente
    //imshow("Imagen Original", img);
    //imshow("Imagen Reconstruida - Tap 1X", reconstructed);
    
    waitKey(0);
    return 0;
}


