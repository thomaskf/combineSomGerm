#include<cstdlib>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<fstream>
#include<iostream>
#include<cstring>

#define VERSION "2.2"

#define FIRST_COL 9
#define EMPTY "./."

// --------------------
// Type of input files
// --------------------

#define MAIN_GERMLINE 0
#define OTHER_GERMLINE 1
#define SOMATIC 2 

using namespace std;

// a single tumor record
class single_record {
public:
    string chr;
    string pos;
    vector<string> record;
    vector<int> type; // 0 - 0/0; 1 - 0/1; 2 - 1/1; 3 - ./. or others
    vector<int> sample;
    string ref;
    string alt;
    string line;
    // the following lines are for germline
    string normal;
    int normal_type; // 0 - 0/0; 1 - 0/1; 2 - 1/1
    void showRecords(ofstream& flog);
    void showRecords(ofstream& flog, int item);
};

// the tumor records from a single file
class file_records {
public:
	int inFileType;
    int numTumors;
    string fileName;
    map<string,single_record> tumor_records;
    file_records(string file_name);
    void loadOtherGermlineFile(bool applyPassFilter);
    void loadSomaticFile(bool applyPassFilter);
    void showRecords();
};

// all tumor records from all files
class all_records {
public:
    vector<string> inFileNames;
    vector<int> inFileTypes;
    vector<int> numTumors;
    string outFileName;
    string outInconsistFileName;
    string outLogFileName;
    vector<file_records*> file_record_list;
    string mainFileName;
    string filterFileList;
    bool onlyKeep00; // only keep normal with 0/0
    bool checkConsist; // to check the records are consistence
    
    // constructor
    all_records();
    
    // destructor
    ~all_records();
    
    // 1. get list of file names
    void getFileNames(int argc, const char** argv);
    
    // 2. load other-germline and all somatic files
    void loadFiles();
    
    // 3. load the main germline result file
    //    and report the overlapping records among all
    void getOverLapRecords();
};

// a single full record
// store the records of all tumor samples
class single_full_germline_record {
public:
    string chr;
    string pos;
    vector<string> record; // the dimension = size of tumor samples
    vector<int> type;
    string ref;
    string alt;
    string line;
    string normal;
    int normal_type; // 0 - 0/0; 1 - 0/1; 2 - 1/1
    
    // check with a file record
    // return true is the record is consistent with this file record
    bool compareFileRecord(file_records* t, map<string,single_record>::iterator& itr, int& mismatch_item, bool checkConsist);

private:
    // compare with another record
    // return true is the record is consistent with this single record
    bool compareSingleRecord(single_record t, int& mismatch_item);
};
