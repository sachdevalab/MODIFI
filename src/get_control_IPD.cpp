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
#include <getopt.h>

using namespace std;


long double *kmer_ipd_sum;
unsigned int *kmer_num_table;
float *kmer_mean_table;
long array_size;

int MAX_KMER_COUNT = 100000;
int up = 7;
int down = 3;

int thread_num;
std::mutex kmer_mutex;

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
};

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
        base[i] = pow(4, k-i-1);
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

struct IPD_Result {
    string refName;
    map<int, float> for_ipd_map;
    map<int, float> rev_ipd_map;
};

struct IPD_Control {
    int kmer_count;
    float control;
};

// load ipd
IPD_Result load_ipd(string raw_ipd){
    map<int, float> for_ipd_map;
    map<int, float> rev_ipd_map;
    ifstream ipd_file;
    ipd_file.open(raw_ipd, ios::in);
    string line;
    string ref;

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
        string base;
        getline(ss, base, ',');
        string coverage;
        getline(ss, coverage, ',');
        string tMean;
        getline(ss, tMean, ',');
        string tErr;
        getline(ss, tErr, ',');
        // cout << refName << "\t" << strand << "\t" << tpl << "\t" << coverage << "\t" << tMean << "\t" << tErr << endl;

        int tpl_int;
        try {
            tpl_int = stoi(tpl);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid tpl argument in " << raw_ipd << " - Full line: " << line << std::endl;
            continue; // Skip this line
        } catch (const std::out_of_range& e) {
            std::cerr << "tpl out of range in " << raw_ipd << " - Full line: " << line << std::endl;
            continue; // Skip this line
        }
        // float tMean_float = stof(tMean);
        float tMean_float;
        try {
            tMean_float = std::stof(tMean);
            // std::cout << "Converted value: " << tMean_float << std::endl;
        } catch (const std::invalid_argument& e) {
            std::cerr << raw_ipd <<" " << tpl <<" Invalid argument: " << tMean << " cannot be converted to float" << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << tMean << " is out of range for float" << std::endl;
        }
        if (strand == "1"){
            for_ipd_map[tpl_int] = tMean_float;
        }
        else{
            rev_ipd_map[tpl_int] = tMean_float;
        }
        ref = refName;

    }
    ipd_file.close();
    IPD_Result result;
    result.for_ipd_map = for_ipd_map;
    result.rev_ipd_map = rev_ipd_map;
    result.refName = ref;
    return result;

}

void update_ipd(string raw_ipd, map<int, IPD_Control> for_control_map, map<int, IPD_Control> rev_control_map, string control_ipd){
    ifstream ipd_file;
    ipd_file.open(raw_ipd, ios::in);
    ofstream control_file;
    control_file.open(control_ipd, ios::out);
    string line;
    string ref;
    float control, ipd_ratio;
    int kmer_count;

    // skip first line
    getline(ipd_file, line);
    control_file << "refName,strand,tpl,base,coverage,tMean,tErr,control,ipd_ratio,kmer_count" << endl;
    while (getline(ipd_file, line)){
        stringstream ss(line);
        // sep line by ,
        string refName;
        getline(ss, refName, ',');
        string strand;
        getline(ss, strand, ',');
        string tpl;
        getline(ss, tpl, ',');
        string base;
        getline(ss, base, ',');
        string coverage;
        getline(ss, coverage, ',');
        string tMean;
        getline(ss, tMean, ',');
        string tErr;
        getline(ss, tErr, ',');
        // cout << refName << "\t" << strand << "\t" << tpl << "\t" << coverage << "\t" << tMean << "\t" << tErr << endl;
        int tpl_int;
        try {
            tpl_int = stoi(tpl);
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid tpl argument in " << raw_ipd << " - Full line: " << line << std::endl;
            continue; // Skip this line
        } catch (const std::out_of_range& e) {
            std::cerr << "tpl out of range in " << raw_ipd << " - Full line: " << line << std::endl;
            continue; // Skip this line
        }
        float tMean_float = stof(tMean);
        if (strand == "1"){
            if (for_control_map.find(tpl_int) != for_control_map.end()){
                control = for_control_map[tpl_int].control;
                kmer_count = for_control_map[tpl_int].kmer_count;
                ipd_ratio = tMean_float/control;
                control_file << refName << "," << strand << "," << tpl << ","<< base << "," << coverage << "," << tMean << "," << tErr << "," << control << "," << ipd_ratio << "," << kmer_count << endl;
            }
            // else{
            //     cout << "No control value for " << tpl_int << endl;
            // }
        }
        else{
            if (rev_control_map.find(tpl_int) != rev_control_map.end()){
                control = rev_control_map[tpl_int].control;
                kmer_count = rev_control_map[tpl_int].kmer_count;
                ipd_ratio = tMean_float/control;
                control_file << refName << "," << strand << "," << tpl << ","<< base << "," << coverage << "," << tMean << "," << tErr << "," << control << "," << ipd_ratio << "," << kmer_count << endl;
            }
        }
        ref = refName;

    }
    ipd_file.close();
    control_file.close();
}



int slide_kmer(string fasta_file, Encode encoder, string raw_ipd){
    IPD_Result ipd_result = load_ipd(raw_ipd);
    Fasta fasta;
    string ref_seq = fasta.get_ref_seq(fasta_file);

    int ref_len;
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index;
    time_t t0 = time(0);

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
            comple_kmer_index += encoder.coder[ref_comple_int[j+encoder.k-1-z]]*encoder.base[z];  
        }
        if (all_valid != false){
            int for_real_pos = j + up;
            if (ipd_result.for_ipd_map.find(for_real_pos) != ipd_result.for_ipd_map.end()) {
                // lock_guard<mutex> guard(kmer_mutex); // Lock the mutex
                if (kmer_num_table[kmer_index] < MAX_KMER_COUNT){
                    kmer_ipd_sum[kmer_index] += ipd_result.for_ipd_map[for_real_pos];
                    kmer_num_table[kmer_index] += 1;
                    // cout<< for_real_pos << " "<<kmer_index<< " ipd "<<kmer_ipd_sum[kmer_index]<< " num " << kmer_num_table[kmer_index]  <<endl;
                }
                else{
                    cout << kmer_index<< " <<<too high kmer count>>>  "<<kmer_num_table[kmer_index]<<endl;
                }
            }
            // to do, for reverse ref
            int rev_real_pos = j + down - 1;
            if (ipd_result.rev_ipd_map.find(rev_real_pos) != ipd_result.rev_ipd_map.end()) {
                // lock_guard<mutex> guard(kmer_mutex); // Lock the mutex
                if (kmer_num_table[comple_kmer_index] < MAX_KMER_COUNT){
                    kmer_ipd_sum[comple_kmer_index] += ipd_result.rev_ipd_map[rev_real_pos];
                    kmer_num_table[comple_kmer_index] += 1;
                    // cout<< rev_real_pos << " "<<comple_kmer_index<< " ipd "<<kmer_ipd_sum[comple_kmer_index]<< " num " << kmer_num_table[comple_kmer_index]  <<endl;
                }
                else{
                    cout << comple_kmer_index<< " <<<too high kmer count>>>  "<<kmer_num_table[comple_kmer_index]<<endl;
                }
            }
        }
        // cout << "kmer_index is "<<kmer_index<<endl;
    }        
    delete [] ref_int;
    delete [] ref_comple_int;
    return 0;
}


int map_control(string fasta_file, Encode encoder, string raw_ipd, string control_ipd){
    IPD_Result ipd_result = load_ipd(raw_ipd);
    Fasta fasta;
    string ref_seq = fasta.get_ref_seq(fasta_file);

    map<int, IPD_Control> for_control_map;
    map<int, IPD_Control> rev_control_map;

    int ref_len, m, e;
    char n;
    float control, ipd_ratio;
    unsigned int kmer_index, comple_kmer_index;
    time_t t0 = time(0);

    cout <<"check control..."<<endl;
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
        if (all_valid != false){
            int for_real_pos = j + up;
            if (ipd_result.for_ipd_map.find(for_real_pos) != ipd_result.for_ipd_map.end()) {
                // lock_guard<mutex> guard(kmer_mutex); // Lock the mutex
                    if (kmer_num_table[kmer_index] > 0){
                        control = kmer_mean_table[kmer_index];
                    }
                    else{
                        control = 1;
                    }
                    IPD_Control control_ipd_obj;
                    control_ipd_obj.kmer_count = kmer_num_table[kmer_index];
                    control_ipd_obj.control = control;
                    for_control_map[for_real_pos] = control_ipd_obj;
                    // for_control_map[for_real_pos] = control;
                    // ipd_ratio = ipd_result.rev_ipd_map[for_real_pos]/control;
            }
            // to do, for reverse ref
            int rev_real_pos = j + down - 1;
            if (ipd_result.rev_ipd_map.find(rev_real_pos) != ipd_result.rev_ipd_map.end()) {
                // lock_guard<mutex> guard(kmer_mutex); // Lock the mutex
                    if (kmer_num_table[comple_kmer_index] > 0){
                        control = kmer_mean_table[comple_kmer_index];
                    }
                    else{
                        control = 1;
                    }
                    IPD_Control control_ipd_obj;
                    control_ipd_obj.kmer_count = kmer_num_table[comple_kmer_index];
                    control_ipd_obj.control = control;
                    rev_control_map[rev_real_pos] = control_ipd_obj;
            }
        }
        // cout << "kmer_index is "<<kmer_index<<endl;
    }        
    delete [] ref_int;
    delete [] ref_comple_int;
    update_ipd(raw_ipd, for_control_map, rev_control_map,control_ipd);
    return 0;
}


string get_base_name(string fasta_file) {
    size_t last_slash = fasta_file.find_last_of('/');
    size_t last_dot = fasta_file.find_last_of('.');
    if (last_slash != string::npos && last_dot != string::npos && last_dot > last_slash) {
        return fasta_file.substr(last_slash + 1, last_dot - last_slash - 1);
    }
    return "";
}

void parallele_each_genome(string genome_list_file, Encode encoder, int start_g, int end_g, string ipd_dir){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);

    string each_genome;
    // float match_rate;
    int genome_index = 0;
    while (getline(list_file, each_genome)){
        if (genome_index >= start_g & genome_index < end_g){
            // extract chr name from each_genome name
            string base_name = get_base_name(each_genome);
            string raw_ipd = ipd_dir + "/" + base_name + ".ipd1.csv";
            slide_kmer(each_genome, encoder, raw_ipd);
        }
        genome_index += 1;
    }

    list_file.close();
}

void parallele_each_genome_control(string genome_list_file, Encode encoder, int start_g, int end_g, string ipd_dir, string control_dir){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);

    string each_genome;
    // float match_rate;
    int genome_index = 0;
    while (getline(list_file, each_genome)){
        if (genome_index >= start_g & genome_index < end_g){
            // extract chr name from each_genome name
            string base_name = get_base_name(each_genome);
            string raw_ipd = ipd_dir + "/" + base_name + ".ipd1.csv";
            string control_ipd = control_dir + "/" + base_name + ".ipd2.csv";
            map_control(each_genome, encoder, raw_ipd, control_ipd);
            
        }
        genome_index += 1;
    }

    list_file.close();
}

struct Genome_count {
    int genome_num;
    int each_thread_g_num;
};

Genome_count assign_parallele(string genome_list_file, int thread_num){
    ifstream list_file;
    list_file.open(genome_list_file, ios::in);
    string each_genome;
    int genome_num = 0;
    while (getline(list_file, each_genome)){
        genome_num += 1;
    }
    list_file.close();
    if (genome_num < thread_num){
        thread_num = genome_num;
    }

    int each_thread_g_num = genome_num/thread_num;
    cout<<"There is "<<genome_num<<" genomes."<<endl;
    cout<<"Each thread process "<<each_thread_g_num<<" genomes"<<endl;

    Genome_count genome_count ;
    genome_count.genome_num = genome_num;
    genome_count.each_thread_g_num = each_thread_g_num;
    // return each_thread_g_num;
    return genome_count;
}

// calculate mean of each kmer
void calculate_mean(){
    for (int i = 0; i < array_size; i++){
        if (kmer_num_table[i] != 0){
            kmer_mean_table[i] = kmer_ipd_sum[i] / kmer_num_table[i];
        }
    }
}

void store_kmer_mean(string kmer_mean_file, string kmer_num_file){
    ofstream index_file;
    index_file.open(kmer_mean_file, ios::out | ios::binary);
    index_file.write((char*)&array_size, sizeof(long));
    index_file.write((char*)kmer_mean_table, sizeof(float)*array_size);
    index_file.close();
    // also store kmer_num_table
    ofstream num_file;
    num_file.open(kmer_num_file, ios::out | ios::binary);
    num_file.write((char*)&array_size, sizeof(long));
    num_file.write((char*)kmer_num_table, sizeof(unsigned int)*array_size);
    num_file.close();

}

void load_kmer_mean(string kmer_mean_file, string kmer_num_file){
    // first check if the file exists
    struct stat buffer;
    if (stat (kmer_mean_file.c_str(), &buffer) != 0){
        cout << "kmer_mean_file does not exist." << endl;
        return;
    }

    ifstream index_file;
    index_file.open(kmer_mean_file, ios::in | ios::binary);
    if (index_file.is_open()){
        index_file.read((char*)&array_size, sizeof(long));
        index_file.read((char*)kmer_mean_table, sizeof(float)*array_size);
        index_file.close();
    }
    ifstream num_file;
    num_file.open(kmer_num_file, ios::in | ios::binary);
    if (num_file.is_open()){
        num_file.read((char*)&array_size, sizeof(long));
        num_file.read((char*)kmer_num_table, sizeof(unsigned int)*array_size);
        num_file.close();
    }
    cout << "kmer_mean_table and kmer_num_table are loaded." << endl;
}

int main(int argc, char *argv[]){
    // string ipd_dir = argv[1];
    // string control_dir = argv[2];
    // string fasta_list = argv[3];
    // int thread_num = stod(argv[4]);

    // string ipd_dir = "/home/shuaiw/borg/bench/break/ipd";
    // string control_dir = "/home/shuaiw/borg/bench/break/test/";
    // string fasta_list = "/home/shuaiw/borg/bench/break/contigs_list.txt";
    // int thread_num = 10;  

    // Default values
    string ipd_dir;
    string control_dir;
    string fasta_list;
    int thread_num = 4;
    string given_kmer_mean_file;
    string given_kmer_num_file;

    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--ipd_dir" && i + 1 < argc) {
            ipd_dir = argv[++i];
        } else if (arg == "--control_dir" && i + 1 < argc) {
            control_dir = argv[++i];
        } else if (arg == "--fasta_list" && i + 1 < argc) {
            fasta_list = argv[++i];
        } else if (arg == "--thread_num" && i + 1 < argc) {
            thread_num = stoi(argv[++i]);
        } else if (arg == "--up" && i + 1 < argc) {
            up = stoi(argv[++i]);
        } else if (arg == "--down" && i + 1 < argc) {
            down = stoi(argv[++i]);
        } else if (arg == "--kmer_mean_file" && i + 1 < argc) {
            given_kmer_mean_file = argv[++i];
        } else if (arg == "--kmer_num_file" && i + 1 < argc) {
            given_kmer_num_file = argv[++i];
        }else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            return 1;
        }
    }

    if (ipd_dir.empty() || control_dir.empty() || fasta_list.empty() || thread_num <= 0) {
        cerr << "Usage: " << argv[0] << " --ipd_dir <dir> --control_dir <dir> --fasta_list <file> --thread_num <num>" << endl;
        return 1;
    }


    int k = up + down;
    int start_g = 0;
    int end_g = 0;
    std::vector<std::thread>threads;

    array_size = pow(4, k);

    kmer_ipd_sum = new long double[array_size];
    kmer_num_table = new unsigned int[array_size];
    memset(kmer_ipd_sum, 0, sizeof(long double)*array_size);
    memset(kmer_num_table, 0, sizeof(unsigned int)*array_size);
    kmer_mean_table = new float[array_size];

    Encode encoder;
    encoder.constructer(k);

    Genome_count genome_count = assign_parallele(fasta_list, thread_num);
    if (genome_count.genome_num < thread_num){
        thread_num = genome_count.genome_num;
        cout << "no. of used threads "<< thread_num << endl;
    }
    bool accept_db = false;
    // if given_kmer_mean_file and given_kmer_num_file are not empty and not none, load them
    if (!given_kmer_mean_file.empty() && !given_kmer_num_file.empty() && given_kmer_mean_file != "None" && given_kmer_num_file != "None"){
        // check if the file exists
        struct stat buffer;
        if (stat (given_kmer_mean_file.c_str(), &buffer) != 0){
            cout << "kmer_mean_file does not exist." << endl;
        }
        else if (stat (given_kmer_num_file.c_str(), &buffer) != 0){
            cout << "kmer_num_file does not exist." << endl;
        }
        else{
            accept_db = true;
        }
        if (accept_db == true){
            cout << "Find control IPD db, load them..." << endl;
            cout << given_kmer_mean_file << endl;
            cout << given_kmer_num_file << endl;
            cout << "Please ensure kmer (dp, down) are the same for --up --down and db." << endl;
            load_kmer_mean(given_kmer_mean_file, given_kmer_num_file);
        }
    }
    if (accept_db == false){
        // load kmer_mean_table and kmer_num_table
        string kmer_mean_file = control_dir + "/control_db.up" + to_string(up) + ".down" + to_string(down) + ".mean.dat";
        string kmer_num_file = control_dir + "/control_db.up" + to_string(up) + ".down" + to_string(down) + ".num.dat";
        // load_kmer_mean(kmer_mean_file, kmer_num_file);
        cout << "start run each genome and count kmers..." << endl;
        for (int i=0; i<thread_num; i++){
            start_g = i*genome_count.each_thread_g_num;
            end_g = (i+1)*genome_count.each_thread_g_num;
            if (i == thread_num - 1){
                end_g = genome_count.genome_num+1;
            }
            cout << "Thread " << i << "\t" << start_g <<"\t" << end_g << endl;
            threads.push_back(thread(parallele_each_genome, fasta_list, encoder, start_g, end_g, ipd_dir));
        }
        for (auto&th : threads)
            th.join();
        threads.clear();
        calculate_mean();
        // store_kmer_mean("/home/shuaiw/borg/test/mean.dat", "/home/shuaiw/borg/test/num.dat");
        store_kmer_mean(kmer_mean_file, kmer_num_file);
    }

    start_g = 0;
    end_g = 0;
    for (int i=0; i<thread_num; i++){
        start_g = i*genome_count.each_thread_g_num;
        end_g = (i+1)*genome_count.each_thread_g_num;
        if (i == thread_num - 1){
            end_g = genome_count.genome_num+1;
        }
        // cout << i << "\t" << start_g <<"\t" << end_g << endl;
        threads.push_back(thread(parallele_each_genome_control, fasta_list, encoder, start_g, end_g, ipd_dir, control_dir));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();
    cout << "All done for control IPD inference." << endl;
    return 0;
}


// ./test /home/shuaiw/borg/new_test11/ipd /home/shuaiw/borg/new_test11/control/ /home/shuaiw/borg/new_test11/contigs_list.txt 2
// ./test --ipd_dir /home/shuaiw/borg/new_test11/ipd --control_dir /home
// /shuaiw/borg/new_test11/control/ --fasta_list /home/shuaiw/borg/new_test11/contigs_list.txt --thread_num 2