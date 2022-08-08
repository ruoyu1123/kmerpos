#include <string>
using namespace std;


int read_kmer(const char *k_path, int t);
// search the poses for certain k-mers; Using 64 bit binary int to store the k-mer 

int build_pos(const char  *fasta_file, string out_path, int k, bool type, int t1, int t2);