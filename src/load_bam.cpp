#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <htslib/sam.h>

using namespace std;

unordered_map<int, vector<int>> ipd_signal;

int read_bam(const string& bam) {
    // open the bam file
    samFile *fp_in = hts_open(bam.c_str(), "r");
    if (fp_in == NULL) {
        cerr << "Error: Cannot open bam file" << endl;
        return 1;
    }

    // read the header
    bam_hdr_t *header = sam_hdr_read(fp_in);
    if (header == NULL) {
        cerr << "Error: Cannot read header" << endl;
        hts_close(fp_in);
        return 1;
    }

    // read the alignment
    bam1_t *aln = bam_init1();
    if (aln == NULL) {
        cerr << "Error: Cannot initialize alignment" << endl;
        bam_hdr_destroy(header);
        hts_close(fp_in);
        return 1;
    }

    // read the alignment one by one
    int read_index = 0;
    while (sam_read1(fp_in, header, aln) >= 0) {
        // Read Name
        std::string read_name = bam_get_qname(aln);
        std::cout << "Read: " << read_name << std::endl;
        // Extract IPD (Inter-Pulse Duration) from the "ip" tag (PacBio-specific)
        uint8_t* ipd_data = bam_aux_get(aln, "ip");
        if (ipd_data) {
            ipd_data++; // Move past the tag type character (usually 'B')
            uint32_t num_ipd_values = bam_auxB_len(ipd_data);  // Get array length
            std::cout << num_ipd_values << std::endl;
            // uint32_t* ipd_array = (uint32_t*)bam_auxB_array(ipd_data);  // Get IPD values
            // std::cout << ipd_array << std::endl;

            // std::cout << "IPD values: ";
            // for (uint32_t i = 0; i < num_ipd_values; i++) {
            //     std::cout << ipd_array[i] << " ";
            // }
            // std::cout << std::endl;
        } else {
            std::cout << "No IPD values found for this read." << std::endl;
        }
    }

    // clean up
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    hts_close(fp_in);

    return 0;
}

int main(){
    string bam = "/home/shuaiw/methylation/data/borg/b_contigs/12.align.bam";
    read_bam(bam);
    return 0;
}
//export PKG_CONFIG_PATH=/home/shuaiw/miniconda3/envs/cpp/lib/pkgconfig/