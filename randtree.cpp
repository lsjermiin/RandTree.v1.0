/////////////////////////////////////////////////////////////////////////////////
// Program name       :randtree.cpp
//
// Version            : 1.0
//
// Author             : Lars S Jermiin
//
// Institutions       : Australian National University
//                      Research School of Biology
//                      Acton, ACT 2601, Australia
//
//                      Univerity College Dublin
//                      School of Biology & Environmental Science
//                      Belfield, Dublin 4, Ireland
//
// Emails             : lars.jermiin [at] anu.edu.au
//                      lars.jermiin [at] ucd.ie
//
// URL                : https://github.com/lsjermiin/randtree
//
// Date begun         : 18 June, 2019
//
// Date modified      : 24 October, 2019
//
// Copyright          : Copyright © 2019 Lars Sommer Jermiin.
//                      All rights reserved.
//
// Responsibility     : The copyright holder takes no legal responsibility for
//                      the correctness of results obtained using this program.
//
// Summary            : RandTree is designed to generate a set of random rooted
//                      or unrooted trees with the tips labelled randomly using
//                      names provided by the user.
//
//                      The output is a file with n Newick formatted trees. The
//                      user specifies all input in the command line.
//
/////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

using namespace std;
using std::cin;

int main(int argc, char** argv) {
    unsigned      seed;
    unsigned long trees(0), tips(0), node(0);
    std::string   rooted;
    std::string   Original_Tree("");
    std::string   Edge("(#,#)");
    std::string   inName, outName;
    std::string   str(""), tmp("");
    std::vector<std::string> taxon;
    std::ifstream infile;
    std::ofstream outfile;
    
    if (argc != 4) {
        std::cerr << "\nRandAlign v1.0 Copyright 2019, Lars Jermiin" << std::endl;
        std::cerr << "\nERROR -- use command: randalign <infile> <r|u> <trees>\n" << std::endl;
        std::cerr << "             infile  Text file with a unique sequence label on each line" << std::endl;
        std::cerr << "             r|u     Rooted (r) or unrooted (u) trees" << std::endl;
        std::cerr << "             trees   Number of random trees to generate" << std::endl;
        std::cerr << std::endl;
        exit(1);
    }
    inName = argv[1];
    rooted = argv[2];
    trees = stoi(argv[3]);
    infile.open(inName.c_str());
    if (!infile) {
        std::cerr << "\nInput file not found" << std::endl;
        exit(1);
    }
    if (toupper(rooted[0]) != 'R' && toupper(rooted[0]) != 'U') {
        std::cerr << "\nIncorrect choice of rooting: [r|u]\n" << std::endl;
        exit(1);
    }
    if (trees < 2 && trees > 2000) {
        std::cerr << "\nIncorrect choice of number of trees [2, ... 2000]\n" << std::endl;
        exit(1);
    }
    while (getline(infile, str)) {
        if (!str.empty()) {
            tmp.clear();
            for (std::string::size_type i = 0; i != str.size(); ++i) {
                if (!isblank(str[i])) {
                    tmp.push_back(str[i]);
                }
            }
            taxon.push_back(tmp);
        }
    }
    infile.close();

    // Prepare output file name
    outName.clear();
    for (std::string::size_type i = 0; i != inName.size() && inName[i] != '.'; ++i) {
        outName += inName[i];
    }
    outName = outName + ".nwk";
    
    // Open input file
    outfile.open(outName.c_str());

    // Random number generator - get seed from system clock
    seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
    // Seed Mersenne-Twister random number generator
    std::mt19937_64 generator(seed);

    // Start loop -- for every tree
    for (unsigned int i = 0; i != trees; ++i) {
        if (toupper(rooted[0]) == 'R') {
            tips = 2;
            Original_Tree = "(#,#);";
        } else {
            tips = 3;
            Original_Tree = "(#,#,#);";
        }
        tmp.clear();
        // Start loop -- for every taxon
        for (std::vector<std::string>::size_type j = tips; j != taxon.size(); ++j) {
            // Generate a uniform distribution of random numbers between 0 and j
            std::uniform_int_distribution<unsigned long> distribution(0,j-1);
            // Pick a random number from this distribution
            node = distribution(generator);
            
            // For every character in the current tree
            unsigned count(0);
            for (unsigned long k = 0; k != Original_Tree.size(); ++k) {
                if (Original_Tree[k] != '#') {
                    tmp.push_back(Original_Tree[k]);
                } else {
                    if (count != node) {
                        tmp.push_back(Original_Tree[k]);
                    } else {
                        tmp = tmp + Edge;
                    }
                    ++count;
                }
            }
            Original_Tree = tmp;
            tmp.clear();
        }
        shuffle(taxon.begin(), taxon.end(), generator);
        tmp.clear();
        
        // Populate unlabelled tree with taxon names
        unsigned count(0);
        for (unsigned long k = 0; k != Original_Tree.size(); ++k) {
            if (Original_Tree[k] != '#') {
                tmp.push_back(Original_Tree[k]);
            } else {
                tmp = tmp + taxon[count];
                ++count;
            }
        }
        // print tree to output file
        outfile << tmp << std::endl;
        tmp.clear();
    }
    outfile.close();
    return 0;
}
