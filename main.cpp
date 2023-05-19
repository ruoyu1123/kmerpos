# include "kmerpos.hpp"
# include <iostream>
# include "clipp.h"
# include <vector>
using namespace clipp; using std::cout; using std::string; using std::vector;


int main(int argc, char *argv[]){
    string kmerfile;        // jellyfish result dump
    string fafile;
    string output = "./";          //output path
    int t = 8;
    int k_size = 21;        //k-mer length

    auto cli = (
        option("-t", "--threads").doc("Number of threads [8]") & value("threads", t),
        option("-k", "--kmer").doc("K-mer length [21]") & value("k-size", k_size),
        option("-o").doc("Output path [./]") & value("outpath", output),
        value("input file", kmerfile).doc("Dump file of jellyfish result"),
        value("fasta", fafile).doc("Query fasta file for k-mer")
    );

    auto fmt = doc_formatting{}
		.first_column(8)                           //left border column for text body
		.doc_column(30)                            //column where parameter docstring starts
		.last_column(100);

    if(!parse(argc, const_cast<char **>(argv), cli)) {
		cout << "Usage:\n" << usage_lines(cli, "kmerpos")
     << "\nOptions:\n" << documentation(cli,fmt) << "\nERROR: Required parameter missing\n";
		// throw "Division by zero condition!";
		exit(0);
	}

    #if 0
	if (help){
		cout << "Usage:\n" << usage_lines(cli, "KAfilter", fmt)
     << "\nOptions:\n" << documentation(cli, fmt) << '\n';
    #endif
	
    read_kmer(kmerfile.c_str(),1);
    build_pos(fafile.c_str(), output,k_size, true, 1, t);
    return 0;

    
}
