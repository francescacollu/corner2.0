#include "utility.hpp"

using namespace std;

int main()
{
    corner::matrix m("../../codice_MPO/OUTPUT/LM.L008.m050.Time001000.dat");
    
    vector<double> v0 = m.col(0);
    vector<double> v1 = m.col(4);

    ofstream newFile("provaOUT.txt");

    for(int k = 0; k < v0.size(); k++)
    {
        cout << v0[k] << '\t' << v1[k] << endl;
        newFile << v0[k] << '\t' << v1[k] << endl;
    }

    // std::string s = "scott>=tiger>=mushroom";
    // std::string delimiter = ">=";
    // std::vector<std::string> tokens;

    // size_t pos = 0;
    // std::string token;
    // while ((pos = s.find(delimiter)) != std::string::npos) 
    // {
    //     token = s.substr(0, pos);
    //     std::cout << token << std::endl;
    //     tokens.push_back(token);
    //     s.erase(0, pos + delimiter.length());
    // }
    // std::cout << s << std::endl;
    // std::cout << endl << tokens[0] << '\t' << tokens[1] << endl;


    return 0;
}