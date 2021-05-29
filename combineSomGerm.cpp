#include "combineSomGerm.h"

// show the help menu
void showHelpMenu(int argc, const char** argv) {
    cout << "Syntax: " << argv[0] << " [options]" << endl << endl;
    cout << "Options:" << endl;
    cout << "   To specify the input files" << endl;
    cout << "      -g [main germline result file]" << endl;
    cout << "      -r [other germline result file]" << endl;
    cout << "      -s [combine somatic result file]" << endl << endl;
    cout << "   To specify the output files" << endl;
    cout << "      -o [output vcf file for the overlapping and consistent results]" << endl;
    cout << "      -e1 [output log file with overlapping but inconsistent results]" << endl;
    cout << "      -e2 [output log file with non-overlapping records]" << endl << endl;
    cout << "   To apply the 'PASS' filter on the files" << endl;
    cout << "      -f [file contains the list of files to apply the filter]" << endl << endl;
    cout << "   Other options:" << endl;
    cout << "      -t : Only keep the records with normal cells which are 0/0" << endl;
    cout << "      -c : Only keep the consistent records" << endl << endl;
    cout << "Remark: 1. If there are more than one germline result files (except the main one)" << endl;
    cout << "           or more than one somatic result files, please use the options '-r' or '-s'" << endl;
    cout << "           more than once." << endl;
    cout << "        2. Options '-g', '-o', and '-e' are compulsory." << endl << endl;
    cout << "Example:" << endl;
    cout << "   $./combineSomGerm -g patient16011.gatk.germline.vcf -r patient16011.Monovar.germline.vcf -s patient16011.gatk.somatic.vcf -s patient16011.somaticsniper.somatic.vcf -s patient16011.varscan.somatic.vcf -o patient16011.overlap.vcf -e1 patient16011.inconsist.log -e2 patient16011.err.log -f applyPassFilterFiles.txt" << endl;
}

void tokenizer(string seq, string separators, vector<string>* result) {
	// split the seq into many parts by "separators"
	// the vector<string> *result cannot be NULL
	result->clear();
	int startpos = (int) seq.find_first_not_of(separators);
	while (startpos != (int) string::npos) {
		int endpos = (int) seq.find_first_of(separators, startpos);
		if (endpos != (int) string::npos) {
			result->push_back(seq.substr(startpos, endpos-startpos));
			startpos = (int) seq.find_first_not_of(separators, endpos);
		} else {
			result->push_back(seq.substr(startpos));
			break;
		}
	}
}

int getType(string record) {
    string s;
    if (record.length() >= 3) {
        s = record.substr(0,3);
        if (s == "0/0")
            return 0;
        else if (s == "0/1")
            return 1;
        else if (s == "1/1")
            return 2;
    }
    return 3;
}

void single_record::showRecords(ofstream& flog) {
    int j;
    for (j=0; j<record.size(); j++) {
        flog << " (" << sample[j] << ")" << record[j];
    }
}

void single_record::showRecords(ofstream& flog, int item) {
    int j;
    for (j=0; j<record.size(); j++) {
        if (sample[j] == item) {
            flog << " (" << sample[j]+1 << ")" << record[j];
            break;
        }
    }
}

file_records::file_records(string file_name) {
    fileName = file_name;
}

// load the combine somatic file
void file_records::loadSomaticFile(bool applyPassFilter) {
    string curr_chr;
    string curr_pos;
    string key_str;
    ifstream fin;
    string aline;
    vector<string> tokens;
    single_record tr;
    int i,k;
    int curr_numTumors;
    inFileType = SOMATIC;
    numTumors = 0;
    
    fin.open(fileName.c_str());
    while (getline(fin,aline)) {
        if (aline.length() > 0 && aline[0] != '#') {
            tokenizer(aline, "\t", &tokens);
            if (tokens.size() > FIRST_COL) {
            	// check the PASS filter
            	if (applyPassFilter && tokens[6] != "PASS")
            		continue;
                // number of tumor samples
                curr_numTumors = ((int)tokens.size() - FIRST_COL + 1) / 2;
                if (numTumors == 0)
                    numTumors = curr_numTumors;
                else if (numTumors != curr_numTumors) {
                    cout << "Error in loading the file: " << fileName << endl;
                    cout << "The number of tumor samples is not " << numTumors << endl;
                    cout << aline << endl;
                    exit(1);
                }
                tr.record.clear();
                tr.type.clear();
                tr.sample.clear();
                k=0;
                for (i=FIRST_COL+1; i<tokens.size(); i+=2) {
                    if (tokens[i].length() >= 3 && tokens[i].substr(0,3) != EMPTY) {
                        if (tokens[i].at(1) == '|') {
                            // phasing information
                            if (tokens[i].substr(0,3) == "1|0") {
                                tokens[i].at(0) = '0';
                                tokens[i].at(2) = '1';
                            }
                            tokens[i].at(1) = '/';
                        }
                        tr.record.push_back(tokens[i]);
                        tr.type.push_back(getType(tokens[i]));
                        tr.sample.push_back(k);
                    }
                    k++;
                }
                if (tr.record.size() > 0) {
                    curr_chr = tokens[0];
                    curr_pos = tokens[1];
                    key_str = curr_chr + " " + curr_pos;
                    tr.chr = curr_chr;
                    tr.pos = curr_pos;
                    tr.line = aline;
                    tr.normal = "";
                    tr.normal_type = 3; // no data
                    tumor_records.insert(pair<string,single_record>(key_str, tr));
                }
            }
        }
    }
    fin.close();
}

// load the germline file
void file_records::loadOtherGermlineFile(bool applyPassFilter) {
    string curr_chr;
    string curr_pos;
    string key_str;
    ifstream fin;
    string aline;
    vector<string> tokens;
    single_record tr;
    int i,k;
    int curr_numTumors;
    inFileType = OTHER_GERMLINE;
    numTumors = 0;
    fin.open(fileName.c_str());
    while (getline(fin,aline)) {
        if (aline.length() > 0 && aline[0] != '#') {
            tokenizer(aline, "\t", &tokens);
            if (tokens.size() > FIRST_COL) {
            	// check the PASS filter
            	if (applyPassFilter && tokens[6] != "PASS")
            		continue;
                // number of tumor samples
                curr_numTumors = (int)tokens.size() - FIRST_COL - 1;
                if (numTumors == 0)
                    numTumors = curr_numTumors;
                else if (numTumors != curr_numTumors) {
                    cout << "Error in loading the file: " << fileName << endl;
                    cout << "The number of tumor samples is not " << numTumors << endl;
                    cout << aline << endl;
                    exit(1);
                }
                tr.record.clear();
                tr.type.clear();
                tr.sample.clear();
                k=0;
                for (i=FIRST_COL+1; i<tokens.size(); i++) {
                    if (tokens[i].length() >= 3 && tokens[i].substr(0,3) != EMPTY) {
                        if (tokens[i].at(1) == '|') {
                            // phasing information
                            if (tokens[i].substr(0,3) == "1|0") {
                                tokens[i].at(0) = '0';
                                tokens[i].at(2) = '1';
                            }
                            tokens[i].at(1) = '/';
                        }
                        tr.record.push_back(tokens[i]);
                        tr.type.push_back(getType(tokens[i]));
                        tr.sample.push_back(k);
                    }
                    k++;
                }
                if (tr.record.size() > 0) {
                    curr_chr = tokens[0];
                    curr_pos = tokens[1];
                    key_str = curr_chr + " " + curr_pos;
                    tr.chr = curr_chr;
                    tr.pos = curr_pos;
                    tr.line = aline;
                    tr.normal = tokens[FIRST_COL]; // normal
                    tr.normal_type = getType(tokens[FIRST_COL]);
                    tumor_records.insert(pair<string,single_record>(key_str, tr));
                }
            }
        }
    }
    fin.close();
}

// load the main germline result file
// and report the overlapping records among all
void all_records::getOverLapRecords() {
    string curr_chr;
    string curr_pos;
    string key_str;
    ifstream fin;
    ofstream fout, flog, finconsist;
    string aline;
    vector<string> tokens;
    single_full_germline_record tr;
    int i,j;
    int numTumors = 0;
    int curr_numTumors;
    bool isConsistent;
    map<string,single_record>::iterator itr;
    int mismatch_item;
    int valid_num;
    bool isNormDiffTumor;
    
    // open file for output
    fout.open(outFileName.c_str());
    flog.open(outLogFileName.c_str());
    finconsist.open(outInconsistFileName.c_str());
    
    fin.open(mainFileName.c_str());
    while (getline(fin,aline)) {
      if (aline.length() > 0) {
          if (aline[0] != '#') {
            // cout << "aline=" << aline << endl << flush;
            tokenizer(aline, "\t", &tokens);
            if (tokens.size() > FIRST_COL) {
                // number of tumor samples
                curr_numTumors = (int)tokens.size() - FIRST_COL - 1;
                if (numTumors == 0)
                    numTumors = curr_numTumors;
                else if (numTumors != curr_numTumors) {
                    cout << "Error in loading the file: " << mainFileName << endl;
                    cout << "The number of tumor samples is not " << numTumors << endl;
                    cout << aline << endl;
                    exit(1);
                }
                curr_chr = tokens[0];
                curr_pos = tokens[1];
                key_str = curr_chr + " " + curr_pos;
                tr.chr = curr_chr;
                tr.pos = curr_pos;
                tr.line = aline;
                isNormDiffTumor = false;
                // normal
                i = FIRST_COL;
                tr.normal = tokens[i];
                tr.normal_type = getType(tokens[i]);
                if (tr.normal_type == 3 || (onlyKeep00 && tr.normal_type > 0)) {
                    continue;
                }
                // other tumor records
                tr.record.clear();
                tr.type.clear();
                valid_num = 0;
                for (i=FIRST_COL+1; i<tokens.size(); i++) {
                    if (tokens[i].length() >= 3 && tokens[i].substr(0,3) != EMPTY) {
                        if (tokens[i].at(1) == '|') {
                            // phasing information
                            if (tokens[i].substr(0,3) == "1|0") {
                                tokens[i].at(0) = '0';
                                tokens[i].at(2) = '1';
                            }
                            tokens[i].at(1) = '/';
                        }
                        tr.type.push_back(getType(tokens[i]));
                        valid_num++;
                    } else {
                        tr.type.push_back(3);
                    }
                    tr.record.push_back(tokens[i]);
                }
                if (valid_num == 0)
                    continue;
                
                // check any tumor is different from normal
                for (i=0; i<tr.type.size(); i++) {
                	if (tr.type[i] != 3 && tr.type[i] != tr.normal_type) {
                		isNormDiffTumor = true;
                		break;
                	}
                }
                if (!isNormDiffTumor)
                	continue;
                
                // compare with all germline and somatic records
                isConsistent = true;
                mismatch_item = 0;
                for (j=0; j<file_record_list.size(); j++) {
                    // cout << "j=" << j << endl << flush;
                    if (!tr.compareFileRecord(file_record_list[j], itr, mismatch_item, checkConsist)) {
                        isConsistent = false;
                        break;
                    }
                }
                if (isConsistent) {
                    fout << aline << endl;
                } else if (j<file_record_list.size() && itr != file_record_list[j]->tumor_records.end()) {
                    finconsist << "[" << mainFileName << "] " << aline << endl;
                    finconsist << "[" << file_record_list[j]->fileName << "] " << itr->second.line << endl;
                    finconsist << "[" << mainFileName << "] (" << mismatch_item+1 << ")" << tr.record[mismatch_item] << endl;
                    finconsist << "[" << file_record_list[j]->fileName << "] ";
                    itr->second.showRecords(finconsist, mismatch_item);
                    finconsist << endl << endl;
                } else if (j<file_record_list.size()) {
                    flog << "[" << mainFileName << "] " << aline << endl;
                    flog << "[" << file_record_list[j]->fileName << "] missing" << endl << endl;
                }
            }
          } else {
            fout << aline << endl;
          }
      }
    }
    fin.close();

    // close file for output
    fout.close();
    flog.close();
    finconsist.close();
}


// print out the tumor records
void file_records::showRecords() {
    int j;
    map<string,single_record>::iterator itr;
    single_record tumr;
    
    for (itr=tumor_records.begin(); itr!=tumor_records.end(); itr++) {
        tumr = itr->second;
        cout << tumr.chr << " " << tumr.pos << endl;
        for (j=0; j<tumr.record.size(); j++) {
            if (j>0)
                cout << " ";
            cout << "(" << tumr.sample[j] << ")" << tumr.record[j];
        }
        cout << endl;
    }
}

// constructor
all_records::all_records() {
	onlyKeep00 = false;
	checkConsist = false;
}

// destructor
all_records::~all_records() {
    int i;
    for (i=0; i<file_record_list.size(); i++) {
        delete file_record_list[i];
    }
}

// get the list of file names from the arguments
void all_records::getFileNames(int argc, const char** argv) {
    int i,k;
    string error_message;
    
    // initialize the fset
    inFileNames.clear();
    inFileTypes.clear();
    outFileName = "";
    outInconsistFileName = "";
    outLogFileName = "";
    mainFileName = "";
    filterFileList = "";
    
    i=1;
    k=0;
    while (i<argc) {
    	if (strcmp(argv[i],"-t")==0) {
    		onlyKeep00 = true;
    		i++;
    	} else if (strcmp(argv[i],"-c")==0) {
    		checkConsist = true;
    		i++;
    	} else {
			if (i+1 < argc) {
				if (strcmp(argv[i],"-g")==0) {
					// -g [main germline result]
					mainFileName = argv[i+1];
					k++;
				} else if (strcmp(argv[i],"-r")==0) {
					// -r [other germline result]
					inFileNames.push_back(argv[i+1]);
					inFileTypes.push_back(OTHER_GERMLINE);
				} else if (strcmp(argv[i],"-s")==0) {
					// -s [somatic result]
					inFileNames.push_back(argv[i+1]);
					inFileTypes.push_back(SOMATIC);
				} else if (strcmp(argv[i],"-o")==0) {
					// -o [output vcf file for the overlapped and consistent results]
					outFileName = argv[i+1];
				} else if (strcmp(argv[i],"-e1")==0) {
					// -e1 [output log file showing the overlapped but inconsistent results]
					outInconsistFileName = argv[i+1];
				} else if (strcmp(argv[i],"-e2")==0) {
					// -e2 [output log file showing the non-overlapping results]
					outLogFileName = argv[i+1];
				} else if (strcmp(argv[i],"-f")==0) {
					// -f [file contains the list of files to apply the filter]
					filterFileList = argv[i+1];
				}
			}
			i+=2;
        }
    }
    
    // checking the inputs
    error_message = "";
    if (i==1)
        error_message = "Error! No valid option is provided";
    else if (k==0)
        error_message = "Error! No main germline result file is provided";
    else if (k>1)
        error_message = "Error! More than one main germline result files are provided";
    else if (outFileName == "")
        error_message = "Error! No output file is provided";
    else if (outInconsistFileName == "")
        error_message = "Error! Option -e1 has to be specified";
    else if (outLogFileName == "")
        error_message = "Error! Option -e2 has to be specified";
    if (error_message != "") {
        showHelpMenu(argc, argv);
        cout << endl;
        cout << error_message << endl;
        exit(1);
    }
    
    // output the state
    if (onlyKeep00) {
    	cout << "Only keep the records with normal cells = 0/0" << endl;
    }
}

// load other-germline and all somatic files
void all_records::loadFiles() {
    int i;
    file_records* curr_record;
    set<string> file_list;
    ifstream fin;
    string aline;
    set<string>::iterator itr;
    bool applyPassFilter;
    
    // first load the file list
    if (filterFileList.length() > 0) {
    	fin.open(filterFileList.c_str());
    	while (getline(fin, aline)) {
    		if (aline.length() > 0) {
    			file_list.insert(aline);
    		}
    	}
    	fin.close();
    }
    
    for (i=0; i<inFileNames.size(); i++) {
        curr_record = new file_records(inFileNames[i]);
        itr = file_list.find(inFileNames[i]);
        applyPassFilter = (itr != file_list.end());
        if (inFileTypes[i] == SOMATIC)
            curr_record->loadSomaticFile(applyPassFilter);
        else // i.e. inFileTypes[i] == OTHER_GERMLINE
            curr_record->loadOtherGermlineFile(applyPassFilter);
        file_record_list.push_back(curr_record);
    }
}

// compare with another record
// assume this is a germline record
// return true is the record is consistent with this record
bool single_full_germline_record::compareSingleRecord(single_record t, int& mismatch_item) {
    int i,k;
    bool isConsistent;
    
    isConsistent = true;
    for (i=0; i<t.record.size(); i++) {
        if (t.type[i] < 3) {
            k = t.sample[i];
            if (type[k] != t.type[i]) {
                isConsistent = false;
                mismatch_item = k;
                break;
            }
        }
    }
    return isConsistent;
}

// check with a file record
// return true is the record is consistent with this file record
bool single_full_germline_record::compareFileRecord(file_records* t, map<string,single_record>::iterator& itr, int& mismatch_item, bool checkConsist) {
    string key;
    
    if (normal_type == 0 || t->inFileType < 2) {
        key = chr + " " + pos;
        itr = t->tumor_records.find(key);
        if (itr == t->tumor_records.end()) {
            // this position does not exist in t
            // cout << "this position does not exist in t" << endl;
            return false;
        }
        if (checkConsist)
        	return compareSingleRecord(itr->second, mismatch_item);
    }
    return true;
}

int main(int argc, const char** argv) {
    
    all_records records;
    
    cout << "CombineSomGerm -- Version " << VERSION << endl;
    cout << "This program reports the overlapped results between different germline and somatic results" << endl << endl;
    
    // show the help menu
    if (argc < 3) {
        showHelpMenu(argc, argv);
        exit(1);
    }
    
    // get list of file names from the arguments
    records.getFileNames(argc, argv);
    
    // load files
    records.loadFiles();
    
    // load the main germline result file
    // and report the overlapping records among all
    // cout << "report the overlapping records among all" << endl << flush;
    records.getOverLapRecords();
    
    cout << "Output files:" << endl;
    cout << "1. " << records.outFileName << " -- Overlapped and consistent vcf records" << endl;
    cout << "2. " << records.outInconsistFileName << " -- Log file with overlapped but inconsistent records" << endl;
    cout << "3. " << records.outLogFileName << " -- Log file with non-overlapping records" << endl << endl;
    cout << "finished!" << endl;
}
