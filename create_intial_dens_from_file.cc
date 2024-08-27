#include <iostream>
#include <ostream>
#include <fstream>

using std::cout;
using std::endl;
using namespace std;

//inputs are input file name, output file name
int main(int argc, char *argv[])
{
    
    const std::string& input_file = std::string(argv[1]);
    const std::string& output_file = std::string(argv[2]);
    
    std::ifstream ogfile;
    std::ofstream newfile;
    ogfile.open(input_file);
    newfile.open(output_file);
    
    ogfile.seekg(EOF); // go to the end of file
    while(true)
    {
        ogfile.unget(); //go back one char
        char in = ogfile.get(); //extract character
        if(in == '\n')
        {
            int i = 0;
            std::string mystring;
            while(ogfile){
                //get data one by one from ogfile and put it into newfile
                std::getline(ogfile, mystring, ',');
                //this takes care of first two elements of file being x_initial and dx_initial
                if(i>1){
                    newfile << std::stod(mystring) << ",";
                }
                i++;
            }
        }
    }
    
    newfile.close();
    ogfile.close();
    
    return 0;
    
}
