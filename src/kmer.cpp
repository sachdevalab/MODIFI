#include <iostream>
#include <fstream>
#include <cmath> 
#include <cstdlib> 
#include <string>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <sys/stat.h> 
#include <pthread.h>
#include <thread>
#include <sstream>
#include <map>
#include <mutex>
#include <queue>
#include <algorithm>
#include <list>

using namespace std;



int thread_num;
std::mutex mtx;  

class Encode{
    public:
    char coder [1000];
    char comple [256];
    int base [32];
    short choose_coder[100];
    int k;
    // short *choose_coder; 

    
    void generate_coder(void);
    void generate_complement(void);
    void generate_base(int k);
    void constructer(int given_k);

};

void Encode::constructer(int given_k){
    k = given_k;
    this->generate_coder();
    this->generate_complement();
    this->generate_base(given_k);
    cout <<"Encoder is established."<<endl;
}

void Encode::generate_coder(void){
    // A:65 97 T:116 84 C:99 67 G: 103 71
    // static char coder [1000];
    for (int j = 0; j < 256; j++){
        coder[j] = 5;
    }
    coder[65] = 0;
    coder[97] = 0;

    coder[84] = 1;
    coder[116] = 1;

    coder[67] = 2;
    coder[99] = 2;

    coder[71] = 3;
    coder[103] = 3; 
}

void Encode::generate_complement(void){
    // static char comple [256];
    for (int j = 0; j < 256; j++){
        comple[j] = 0;
    }
    comple[65] = 84;
    comple[97] = 84;
    comple[116] = 65;
    comple[84] = 65;
    comple[99] = 71;
    comple[67] = 71;
    comple[103] = 67;
    comple[71] = 67;
    // return comple;   
}

void Encode::generate_base(int k){
    
    for (int i = 0; i<k; i++){       
        base[i] = pow(2, k-i-1);
    }
    // return base;    
}

string get_read_ID(string reads_seq){
    string delimiter = "/";
    string read_name_forward = reads_seq.substr(0, reads_seq.find(delimiter));
    delimiter = " ";
    read_name_forward = read_name_forward.substr(0, read_name_forward.find(delimiter));
    delimiter = "\t";
    read_name_forward = read_name_forward.substr(0, read_name_forward.find(delimiter));
    return read_name_forward;
}

class Fasta{
    public:
    string ref_seq;
    string line_seq, chr_name;
    string get_ref_seq(string fasta_file);
    string get_ref_name(string fasta_file);
};

string Fasta::get_ref_seq(string fasta_file){
    ifstream fa_file;
    fa_file.open(fasta_file, ios::in);
    // string line_seq;
    ref_seq = "\0";
    while (getline(fa_file,line_seq)){
        if (line_seq[0] == '>'){
            continue;
        }
        else{
            ref_seq += line_seq;           
        }            
    }
    return ref_seq;
}

string Fasta::get_ref_name(string fasta_file){
    ifstream fa_file;
    fa_file.open(fasta_file, ios::in);
    
    while (getline(fa_file,line_seq)){
        if (line_seq[0] == '>'){
            chr_name = get_read_ID(line_seq).substr(1);   
            break; 
        }
    }
    return chr_name;
}

int get_fitness(string fasta_file, Encode encoder){
    Fasta fasta;
    string ref_seq = fasta.get_ref_seq(fasta_file);
    double match_base_num = 0;

    int ref_len;
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index, real_index;
    time_t t0 = time(0);
    int covert_num, comple_num;
    short convert_ref[300];

    // cout <<"check fitness..."<<endl;
    ref_len= ref_seq.length();
    // cout <<"genome len is "<<ref_len<<endl;
 
    int *ref_int = new int[ref_len];
    int *ref_comple_int = new int[ref_len];
    for (int j = 0; j < ref_len; j++){
        ref_int[j] = (int)ref_seq[j];
        ref_comple_int[j] = encoder.comple[ref_int[j]];
    }
    for (int j = 0; j < ref_len-encoder.k+1; j++){
        int max_dp = 0;
        kmer_index = 0;
        comple_kmer_index = 0;
        bool all_valid = true;
        for (int z = 0; z < encoder.k; z++){
            m = encoder.coder[ref_int[j+z]];
            if (m == 5){
                all_valid = false;
                break;
            }
            kmer_index += m*encoder.base[z]; 
            comple_kmer_index += encoder.coder[ref_comple_int[j+z]]*encoder.base[(encoder.k-1-z)];  
        }
        if (all_valid == false){
            kmer_index = -1;
            comple_kmer_index = -1;
        }
        cout << "kmer_index is "<<kmer_index<<endl;
    }        
    delete [] ref_int;
    delete [] ref_comple_int;
    return 0;
}

// load ipd
void load_ipd(string raw_ipd){
    map<int, vector<int>> for_ipd_map;
    map<int, vector<int>> rev_ipd_map;
    ifstream ipd_file;
    ipd_file.open(raw_ipd, ios::in);
    string line;

    // skip first line
    getline(ipd_file, line);

    while (getline(ipd_file, line)){
        stringstream ss(line);
        // sep line by ,
        string refName;
        getline(ss, refName, ',');
        string strand;
        getline(ss, strand, ',');
        string tpl;
        getline(ss, tpl, ',');
        string coverage;
        getline(ss, coverage, ',');
        string tMean;
        getline(ss, tMean, ',');
        string tErr;
        getline(ss, tErr, ',');
        // cout << refName << "\t" << strand << "\t" << tpl << "\t" << coverage << "\t" << tMean << "\t" << tErr << endl;

        int tpl_int = stoi(tpl);
        int tMean_int = stoi(tMean);
        if (strand == "0"){
            for_ipd_map[tpl_int].push_back(tMean_int);
        }
        else{
            rev_ipd_map[tpl_int].push_back(tMean_int);
        }

    }
    ipd_file.close();
}

int main(){
    string fasta_file = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/12.fa";
    string raw_ipd = "/home/shuaiw/methylation/data/borg/b_contigs/test/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_219069_438138.ipd1.csv";

    int up = 7;
    int down = 3;
    int k = up + down;
    Encode encoder;
    encoder.constructer(k);
    // get_fitness(fasta_file, encoder);
    load_ipd(raw_ipd);
    return 0;
}