#include <iostream>
#include <fstream>  // Required for std::ofstream
#include <string>
#include <vector>
#include <unordered_map>
#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordView.h>

using namespace std;
using namespace PacBio::BAM;

// unordered_map<int64_t, vector<uint16_t>> ipd_by_locus;

unordered_map<int64_t, vector<uint16_t>> process_read(const BamRecord& record, unordered_map<int64_t, vector<uint16_t>>& ipd_by_locus) {
    int32_t ref_pos = record.ReferenceStart();
    int32_t query_pos = 0;

    const auto& cigar_data = record.CigarData();
    for (const auto& cigar : cigar_data) {
        PacBio::Data::CigarOperationType cigar_op = cigar.Type();
        uint32_t cigar_len = cigar.Length();

        switch (cigar_op) {
            case PacBio::Data::CigarOperationType::ALIGNMENT_MATCH: // alignment match (can be a sequence match or mismatch)
            case PacBio::Data::CigarOperationType::SEQUENCE_MATCH: // sequence match
            case PacBio::Data::CigarOperationType::SEQUENCE_MISMATCH: // sequence mismatch
                for (uint32_t i = 0; i < cigar_len; ++i) {
                    // cout << "<<<" << ref_pos << "\t" << query_pos << "\t" << record.IPD().Data().size() << endl;

                    if (query_pos < record.IPD().Data().size() && record.IPD().Data()[query_pos] != 0) {
                        ipd_by_locus[ref_pos].push_back(record.IPD().Data()[query_pos]);
                    }
                    // cout << ref_pos << "\t" << query_pos << "\t"  << endl;
                    ++ref_pos;
                    ++query_pos;
                }
                break;
            case PacBio::Data::CigarOperationType::INSERTION: // insertion to the reference
                query_pos += cigar_len;
                break;
            case PacBio::Data::CigarOperationType::DELETION: // deletion from the reference
                ref_pos += cigar_len;
                break;
            case PacBio::Data::CigarOperationType::REFERENCE_SKIP: // skipped region from the reference
                ref_pos += cigar_len;
                break;
            case PacBio::Data::CigarOperationType::SOFT_CLIP: // soft clipping (clipped sequences present in SEQ)
                query_pos += cigar_len;
                break;
            case PacBio::Data::CigarOperationType::HARD_CLIP: // hard clipping (clipped sequences NOT present in SEQ)
                break;
            case PacBio::Data::CigarOperationType::PADDING: // padding (silent deletion from padded reference)
                break;
            default:
                break;
        }
    }
    return ipd_by_locus;
}

int avergage_ipd(unordered_map<int64_t, vector<uint16_t>>& ipd_by_locus, const string& raw_ipd_file) {
    // write raw ipd to file
    ofstream raw_ipd_stream(raw_ipd_file);
    
    // for each locus, calculate the mean and std of IPD
    for (auto& locus : ipd_by_locus) {
        int64_t locus_id = locus.first;
        vector<uint16_t> ipd_values = locus.second;
        int ipd_sum = 0;
        for (auto& ipd : ipd_values) {
            ipd_sum += ipd;
        }
        double ipd_mean = ipd_sum / ipd_values.size();
        double ipd_std = 0;
        for (auto& ipd : ipd_values) {
            ipd_std += (ipd - ipd_mean) * (ipd - ipd_mean);
        }
        ipd_std = sqrt(ipd_std / ipd_values.size());
        // cout << locus_id << "\t" << ipd_mean << "\t" << ipd_std << "\t" << ipd_values.size() << endl;
        raw_ipd_stream << locus_id << "\t" << ipd_mean << "\t" << ipd_std << "\t" << ipd_values.size() << endl;
    }
    raw_ipd_stream.close();
    return 0;
}         

int read_bam(const string& bam, const string& fasta, const string& raw_ipd_file) {
        BamFile bamFile(bam);
        BamReader reader(bamFile);
        BamRecord record;
        unordered_map<int64_t, vector<uint16_t>> ipd_by_locus;

        int read_index = 0;
        while (reader.GetNext(record)) {
            if (! record.IsMapped()){
                continue;
            }
            // Read Name
            std::string read_name = record.FullName();
            std::cout << "Read: " << read_name << std::endl;

            // cal read length
            int read_length = record.Sequence().size();

            // Get the read strand
            PacBio::Data::Strand strand = record.AlignedStrand();
            if (strand == PacBio::Data::Strand::FORWARD) {
                std::cout << "Strand: FORWARD" << std::endl;
                continue;
            } else {
                std::cout << "Strand: REVERSE" << std::endl;
            }
            
            std::cout << "Strand: " << strand << std::endl;

            // Get the alignment locus on the reference
            // int32_t ref_start = record.ReferenceStart();
            // int32_t ref_end = record.ReferenceEnd();
            // std::cout << "Alignment locus on reference: " << ref_start << "-" << ref_end << std::endl;

            // // Extract IPD (Inter-Pulse Duration) from the record
            if (record.HasIPD()) {
                const auto& ipd_values = record.IPD().Data();
                if (ipd_values.size() > 0) {
                    ipd_by_locus = process_read(record, ipd_by_locus);
                }
                


            }
            read_index += 1;
            // if (read_index > 10) {
            //     break;
            // }
            std::cout << "read index" << read_index << "\t" << ipd_by_locus.size() << std::endl;
        } 
        // for each locus, calculate the mean and std of IPD
        int x = avergage_ipd(ipd_by_locus, raw_ipd_file);

    return 0;
}



int main() {
    string bam = "/home/shuaiw/methylation/data/borg/b_contigs/12.align.bam";
    string fasta_file = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/12.fa";
    string raw_ipd_file = "/home/shuaiw/methylation/data/borg/b_contigs/cpp/12_ipd.txt";
    int result = read_bam(bam, fasta_file, raw_ipd_file);
    if (result != 0) {
        cerr << "Error reading BAM file" << endl;
        return 1;
    }


    return 0;
}