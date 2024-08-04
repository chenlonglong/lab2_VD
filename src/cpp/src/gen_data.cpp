#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <queue>
#include <stdlib.h>
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
	uint32_t 	read_len 	= 64;
    uint32_t 	hap_len 	= 128;
	uint32_t 	hap_amount	= 100;
	uint32_t 	read_amount	= 100;
	uint32_t 	data_amount = 100;

    for(int i=0;i<argc;i++) {
        std::string argv_str(argv[i]);

		// Turn the argv into uint32_t (stoul) and int (stoi)		
		if(argv_str == "-read_len") read_len = std::stoul(argv[i+1]);
		if(argv_str == "-hap_len") 	hap_len = std::stoul(argv[i+1]);
		// if(argv_str == "-hap_amount") hap_amount = std::stoul(argv[i+1]);
		// if(argv_str == "-read_amount") read_amount = std::stoul(argv[i+1]);
		if(argv_str == "-data_amount") data_amount = std::stoul(argv[i+1]);
    }
	
	std::cerr << "\n";
	std::cerr << "========= Start Generating dataset =========" << "\n";
    std::cerr << "Read Length: " 		<< read_len 	<< "\n";
	std::cerr << "Hap Length: " 		<< hap_len 		<< "\n";
	std::cerr << "Haplotype amount: " 	<< hap_amount 	<< "\n";
	std::cerr << "Read amount: " 		<< read_amount 	<< "\n";
	std::cerr << "data_amount: "        << data_amount	<< "\n";
	std::cerr << "========================================" << "\n\n";

	if (read_len >= hap_len) {
		std::cerr << "Read length has to be smaller than haplotype length! \n";
		return 1;
	}
	if ((read_len % 4) != 0) {
		std::cerr << "Read length % 4 should be 0! \n";
		return 1;
	}
	if ((hap_len % 4) != 0) {
		std::cerr << "Haplotype length % 4 should be 0! \n";
		return 1;
	}
	if (data_amount < 100) {
		std::cerr << "data amount should be greater than 100!\n";
	}

	std::ofstream fout_data("./test_data/random_pattern.txt");
	std::ofstream fout_bin("./test_data/random_pattern.bin", ios::out | ios::binary);
	if(!fout_bin) {
      cout << "Cannot open binary file!" << endl;
      return 1;
   	}

	srand (2);		//random seed
	int random_num;	// random number
	std::string haplotype_seq, read_seq, replace_seq;
	std::vector<string> haplotypes, reads;

	char* hap_buffer  = new char[32];
	char* read_buffer = new char[32];
	char* read_BQ_buffer = new char[32];
	int temp_byte = 0;

	std::string reference;

	for(int i = 0; i < 1000000; i++) {
		random_num = rand();

		if ((random_num % 4) == 0)		reference.append("A");
		else if ((random_num % 4) == 1)	reference.append("C");
		else if ((random_num % 4) == 2)	reference.append("G");
		else if ((random_num % 4) == 3)	reference.append("T");
	}
	
	for (int i = 0; i < hap_amount; i++) {
		auto hap_offset = rand() % (1000000 - hap_len);
		haplotype_seq = reference.substr(hap_offset, hap_len);
		
		random_num = rand();
		// if (random_num % 4 == 0) // exact match do nothing
		if (random_num % 4 == 1) {
			auto mismatch_amount = rand() % 10 + 1;
			for (int m = 0; m < mismatch_amount; m++) {
				int offset = rand() % haplotype_seq.size();
				
				if (haplotype_seq[offset] == 'A')	haplotype_seq.replace(offset, 1, "C");
				else								haplotype_seq.replace(offset, 1, "A");
			}
		}
		else if (random_num % 4 == 2) { // deletion
			auto delete_len = rand() % 20 + 1;
			auto delete_pos = rand() % (hap_len - delete_len);
			haplotype_seq.erase(delete_pos, delete_len);
			haplotype_seq = haplotype_seq + reference.substr(hap_offset + hap_len, delete_len);
		}
		else if (random_num % 4 == 3) { // insertion
			auto insert_len = rand() % 20 + 1;
			auto insert_pos = rand() % (hap_len - insert_len);
			for (int k = insert_pos; k < insert_len; k++) {
				haplotype_seq.insert(insert_pos, "A");
			}
		}



		haplotypes.push_back(haplotype_seq);
		haplotype_seq = "";
	}

	for (int i = 0; i < read_amount; i++) {
		int hap_id = rand() % hap_amount;
		int read_offset = rand() % (hap_len - read_len);
		
		read_seq =  haplotypes[hap_id].substr(read_offset, read_len);
		reads.push_back(read_seq);
	}


	// For binary
	for (int i = 0; i < data_amount; i++){
		int hap_idx = rand() % hap_amount;
		int read_idx = rand() % read_amount;
		// Write binary file (Ref)
		auto hap_temp = haplotypes[hap_idx];
		auto read_temp = reads[read_idx];
		std::string read_BQ_temp;
		for (int j = 0; j < read_len; j++) read_BQ_temp.append("3");
		std::cout << "hap    : " << hap_temp << std::endl;
		std::cout << "read   : " << read_temp << std::endl;
		std::cout << "read BQ: " << read_BQ_temp << std::endl << std::endl;

		for (int k = 0; k < hap_len; k++){
			if      (hap_temp[k] == 'A')	temp_byte += 0;
			else if (hap_temp[k] == 'C')	temp_byte += 1;
			else if (hap_temp[k] == 'G')	temp_byte += 2;
			else if (hap_temp[k] == 'T')	temp_byte += 3;

			if (k % 4 == 3) {
				hap_buffer[k/4] = temp_byte;
				temp_byte = 0;
			}
			else temp_byte = temp_byte * 4;
		}
		fout_bin.write (hap_buffer, (hap_len/4));

		// Write binary file (Read)
		temp_byte = 0;
		for (int k = 0; k < read_len; k++){
			if      (read_temp[k] == 'A')  temp_byte += 0;
			else if (read_temp[k] == 'C')	temp_byte += 1;
			else if (read_temp[k] == 'G')	temp_byte += 2;
			else if (read_temp[k] == 'T')	temp_byte += 3;

			if (k % 4 == 3) {
				read_buffer[k/4] = temp_byte;
				temp_byte = 0;
			}
			else temp_byte = temp_byte * 4;
		}
		fout_bin.write (read_buffer, (read_len/4));

		// Write binary file (Read)
		temp_byte = 0;
		for (int k = 0; k < read_len; k++){
			if      (read_BQ_temp[k] == '0')  temp_byte += 0;
			else if (read_BQ_temp[k] == '1')	temp_byte += 1;
			else if (read_BQ_temp[k] == '2')	temp_byte += 2;
			else if (read_BQ_temp[k] == '3')	temp_byte += 3;

			if (k % 4 == 3) {
				read_BQ_buffer[k/4] = temp_byte;
				temp_byte = 0;
			}
			else temp_byte = temp_byte * 4;
		}
		fout_bin.write (read_BQ_buffer, (read_len/4));






		fout_data << hap_temp  << "\n";
		fout_data << read_temp << "\n";
		fout_data << read_BQ_temp << "\n\n";
	}





	
	fout_data.close();
	fout_bin.close();

	delete[] read_buffer;
	delete[] hap_buffer;

    return 0;
}
