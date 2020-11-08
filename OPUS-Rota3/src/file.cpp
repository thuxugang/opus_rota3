#include <fstream>
#include <iostream>
#include "file.h"
#include "config.h"

using namespace std;


vector<string> getFiles(string& path){
    vector<string> files;
    ifstream inputList;
    inputList.open(path.c_str());
    string str_line;
    if (inputList.is_open()){
        while (!inputList.eof()){
            getline(inputList, str_line);
            if (str_line == "" || str_line == " " || str_line == "\t" || str_line == "\r"){
                continue;
            }
            files.push_back(trim(str_line));
        }
    }else{
        cout << "Cannot open input list" << endl;
        exit(-1);
    }
	return files;
}

