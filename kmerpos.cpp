#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <bitset>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include "thread_pool.hpp"
#include "uthash.h"
#include "kmerpos.hpp"
#include "stringtools.hpp"

using namespace std;

// vector<phmap::flat_hash_set<uint64_t>> kmer; //存储kmer的hash>
mutex m;
int k_size;          // k-mer len
unordered_set<uint64_t> kmer;

#if 0
uthash table
this hashtable need more memeories to store pair;

typedef struct my_struct {
    uint64_t key;
    //char name;
    UT_hash_handle hh;
} hash_node;

class hash_set{
public:
    hash_node *dict=NULL;
    void insert(uint64_t *key){
        hash_node *s;
        HASH_FIND_INT(dict, key, s);
        if (s == NULL)
        {
            s = (struct my_struct*)malloc(sizeof *s);
            s->key = *key;
            HASH_ADD_INT(dict, key, s);
        }
    }
    hash_node *find(uint64_t *key)
    {
        hash_node *tmp=NULL;
        HASH_FIND_INT(dict, key, tmp);
        return tmp;
    }
    int size() {
        return HASH_COUNT(dict);
    }
};
vector<hash_set> kmer;
#endif


std::unordered_map<char, uint64_t> ctoi = {
	{'A', 0ull},
	{'T', 3ull},
	{'C', 1ull},
	{'G', 2ull}
};

/*
// 读取k-mer为二进制整数
uint64_t ctoi(char c)
{
    switch (c)
    {
    case 'A':
        return 0ull;
    case 'T':
        return 3ull;
    case 'C':
        return 1ull;
    case 'G':
        return 2ull;
    default:
        return 4ull;
    }
}
*/



//
uint64_t tobin(string* s){
    uint64_t k_value = 0;
    int i = 0;
    for (auto c = s->begin(); c < s->end(); c++){
        uint64_t v = ctoi[*c];
        k_value = k_value << 2 | v;
    }
    return k_value;
}


// reverse the bin of kmer
//uint64_t reversebin(uint64_t x){
//	return;
//}

int read_kmer(const char *k_path, int k)
{
    ifstream k_mer(k_path);
    uint64_t i = 0;
    uint64_t v;
    uint64_t k_value = 0;
    // printf("%s",s);
	string line;
    while(std::getline(k_mer, line))
    {
		string tmp = line.substr(0, k);	
		v = tobin(&tmp);
		kmer.insert(v);
	}
	printf("Total kmer number: %ld", kmer.size());
    return 1;
}

/* store the sbin(bin of kmer),rbin(reverse bin of kmer), pos at contig and contig name*/
class KMER
{
public:
    uint64_t sbin, rbin;
    uint32_t pos;
    int nameid;
    int is_ukmer(FILE *fp)
    {
        if (kmer.find(sbin) != kmer.end())
        {
            // cout<<nameid<<"\t"<<pos<<"\t"<<sbin<<"\t+"<<endl;
            fprintf(fp, "%x\t%lx\n", pos, sbin);
            // return '+';
        }
        else if (kmer.find(rbin) != kmer.end())
        {
            // cout<<nameid<<"\t"<<"\t"<<i<<"\t"<<rbin<<"\t-\n";
            fprintf(fp, "%x\t%lx\n", pos, rbin);
            // return '-';
        }
        return 0;
    }
};

int search_kmer(string line, int t1, string name, FILE *file[], bool mask[], int n)
{
    // Find a idle file
    FILE *fp;
    int index;
    m.lock();
    for (int i = 0; i < n; i++)
    {
        if (mask[i] == 0)
        {
            index = i;
            fp = file[i];
            mask[i] = 1;
            break;
        }
    }
    m.unlock();
    KMER k;
    vector<string> id = split_find(name, " ");
    fprintf(fp, "@%s\n", id[0].c_str());
    uint64_t base;
    uint32_t len = line.length();
    string qkmer = line.substr(0, k_size);
    k.sbin = tobin(&qkmer);
    reverse(qkmer.begin(),qkmer.end());
    k.rbin = ~tobin(&qkmer);
    cout<<k.rbin<<" "<<k.sbin<<"\n";
    for (uint32_t i = k_size; i < len; i++)
    {
        k.pos = i - (k_size-1);
        k.is_ukmer(fp);
        base = ctoi[line[i]];
        k.sbin = k.sbin << 2 | base & 0x3ffffffffffull;
        k.rbin = k.rbin >> 2 | ~base << 40;
    }
    k.is_ukmer(fp);
    m.lock();
    mask[index] = 0;
    m.unlock();
    return 1;
}

int build_pos(const char *fasta_file, string out_path,int k, bool type, int t1, int t2)
{
    cout << "Number of threads: " << t2 << "\n"
         << "Reading the kmer to dicts\n";
    ThreadPool fpool(t2);
    fpool.init();
    ifstream fa(fasta_file);
    string line, name;
    k_size = k;
    FILE *file[t2];

    /*
        Creat mutiple files for mutiple threads;
        File pointers are store in a FILE * array.
    */
    char outname[15];
    if (type)
        strcpy(outname, "/kmerAA_r.pos");
    else
        strcpy(outname, "/kmerAA_q.pos");
    string mk = "mkdir -p " + out_path;
    system(mk.c_str());

    for (int i = 0; i < t2; i++)
    {
        char *outfile = (char *)malloc((out_path.size() + 25) * sizeof(char));
        strcpy(outfile, out_path.c_str());
        strcat(outfile, outname);
        file[i] = fopen(outfile, "w");
        if (outname[5] == 'Z')
            outname[5] = 'a';
        else
            outname[5] = outname[5] + 1;
        if (outname[5] == 'z')
        {
            if (outname[6] == 'Z')
                outname[6] = 'a';
            else
                outname[6]++;
        }
    }
    bool mask[t2] = {0};            // In order to search a idle file to write
    // auto size=sizeof(KMER);
    cout << "Searching kmer in fasta file!\n";
    while (getline(fa, line))
    {
        if (line[0] == '>')
        {
            name = line.substr(1, line.find(" ")-1);
        }
        else
            fpool.submit(search_kmer, line, t1, name, &file[0], &mask[0], t2);
    }
    fpool.shutdown();
    for (int i = 0; i < t2; i++)
    {
        fclose(file[i]);
    }
    fa.close();
    string cmd = "cat " + out_path + "/*_r.pos > " + out_path + "/ref.pos;" + "rm " + out_path + "/*_r.pos";
    system(cmd.c_str());
#if 0
    string cmd,cmd2;
	if(type){
		cmd="cat " + out_path +"/*_r.pos > "+ out_path + "/ref.pos";
		cmd2="rm "+out_path + "/*_r.pos";
	}
	else{
		cmd="cat " + out_path +"/*_q.pos > "+ out_path + "/query.pos";
		cmd2="rm "+out_path + "/*_q.pos";
	}
	system(cmd.c_str());
	system(cmd2.c_str());
#endif

    return 1;
}
