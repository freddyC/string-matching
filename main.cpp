// Fred Christensen

#include "tools.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <random>

typedef std::map< std::string, std::map <std::string, int> > nextWordCount;
typedef std::map< std::string, std::pair <std::vector< std::string >, std::vector< int > > > nextWordProb;

////////////// PATTERN MATCHING ALGORITHMS //////////////
std::vector<int> nieve (std::string, std::string);

std::vector<int> calcJumpTable (std::string pat);
std::vector<int> kml (std::string, std::string, std::vector<int>);

std::map <char, int> calcDeltaTable (std::string);
std::vector<int> bm (std::string, std::string, std::map<char, int>);

////////////// STRING GENERATING (MARCOV'S CHAIN) //////////////
nextWordProb calcFile (std::string);
void generateText (std::string, nextWordProb, unsigned long long int);

////////////  HELPER FUNCTIONS //////////////
nextWordCount mapWord (nextWordCount, std::string, std::string);
std::string randomKey (nextWordProb);
std::string findNext (nextWordProb, std::string);
nextWordProb calcNextProb (nextWordCount nwc);
std::string readFile (std::string filename);
void testAlgorithms (std::string txt, std::string pat);
void printRes (std::vector<int>);
void printStr (std::string res);

int main() {

  std::string inFile = "/Users/fred.christensen/Dropbox/school/Algorithms/string-matching/dna.txt";
  std::string outFile = "/Users/fred.christensen/Dropbox/school/Algorithms/string-matching/new-dna.txt";

  auto nwp = calcFile(inFile);
  generateText(outFile, nwp, 4);

  std::cout << "String Matching Comparisons\n";
  std::string txt = readFile(inFile);
  std::string pat = readFile(outFile);
  testAlgorithms(txt, pat);

  return 0;
}

std::string readFile (std::string filename) {
  std::ifstream fin (filename.c_str());
  std::string str;
  std::string temp;
  while (fin >> temp) {
    str.append(temp);
  }
  return str;
}

void testAlgorithms(std::string txt, std::string pat) {
  std::vector<int> nieveRes = nieve(txt, pat);

  std::map <char, int> delta = calcDeltaTable(pat);
  std::vector<int> bmRes = bm(txt, pat, delta);

  std::vector<int> jumpTable = calcJumpTable(pat);
  std::vector<int> kmlRes = kml(txt, pat, jumpTable);

  bool areEqual = true;
  if (bmRes.size() != kmlRes.size() || kmlRes.size() != nieveRes.size()) {
    areEqual = false;
  }
  for (int i = 0; i < bmRes.size() && areEqual; ++i) {
    if (bmRes[i] != kmlRes[i] || kmlRes[i] != nieveRes[i]) {
      areEqual = false;
    }
  }
  auto res = (areEqual) ? "TRUE\n\n" : "FALSE\n\n";
  int length = std::max( nieveRes.size(), bmRes.size());
  length = std::max( nieveRes.size(), kmlRes.size());
  std::cout << "result size:\t" << length;
  std::cout << "\nresult consistent:\t" << res;
}

void printRes (std::vector<int> res) {
  for (int i : res) {
    std::cout << i << ", ";
  }
  std::cout << std::endl;
}

void printStr (std::string res) {
  for (char i : res) {
    std::cout << i << ", ";
  }
  std::cout << std::endl;
}

std::vector<int> nieve (std::string text, std::string pat) {
  std::vector<int> results;
  for(int i = pat.size() -1; i < text.size(); ++i) {
    bool matches = true;
    for (int j = pat.size() -1, k = 0; j >= 0; --j, ++k) {
      if (text[i - k] != pat [j]) {
        j = -1;
        matches = false;
      }
    }
    if (matches == true) {
      results.push_back(i - pat.size() +1);
    }
  }
  return results;
}

// Helped by http://ankitsjain22.wordpress.com/2013/07/17/kmp-algorithm-c-code/
std::vector<int> calcJumpTable (std::string pat) {
  std::vector<int> jumpTable (pat.size());

  // indexes 0 and 1 will always be 0
  jumpTable[0] = 0;
  jumpTable[1] = 0;

  int cur = 0;
  for (int i = 2; i < pat.size(); ++i) {
    while(cur != 0 && pat[cur] != pat[i-1]) {
      cur = jumpTable[cur];
    }

    if(pat[cur] == pat[i-1]) {
      cur++;
    }

    jumpTable[i] = cur;
  }
  return jumpTable;
}

std::vector<int> kml (std::string text, std::string pat, std::vector<int> jumpTable) {
  std::vector<int> res;
  int cur = 0;
  for(int i = 0; i<text.size(); i++) {
    while(cur > 0 && pat[cur] != text[i]) {
      cur = jumpTable[cur];
    }

    if(pat[cur] == text[i]){
      if(++cur == pat.size()) {
        res.push_back(i - pat.size() +1);
        cur = 0;
      }
    }
  }

  return res;
}


// this is all me and it is beutiful :)
std::map <char, int> calcDeltaTable (std::string pat) {
  std::map <char, int> delta;
  for (int i = 0; i < pat.size(); ++i) {
    delta[pat[i]] = pat.size() -i;
  }
  return delta;
}

// Followed patterns found here
// http://www.geeksforgeeks.org/pattern-searching-set-7-boyer-moore-algorithm-bad-character-heuristic/
std::vector<int> bm (std::string text, std::string pat, std::map<char, int> delta) {
  std::vector<int> results;
  for (int i = 0; i <= text.size() - pat.size(); ++i) {
    int j = pat.size() -1;
    for (; j >= 0 && pat[j] == text[i+j];--j){
      if (j == 0) {
        results.push_back(i);
      }
    }
    if (j < 0) {
      i += std::max(1, j - delta[text[i+j]]);
    }
  }
  return results;
}


// CREATING BIAS FOR RANDOMIZATION
// read file
// map each word to a vector of pairs
//    each pair has a possible following word, and the number of times it follows
// Sudo Code
//    read first word
//    while not at end of file
//      read next word
//      map first word to next word
//      first word = next word
//    calc probability map
//    return probability map

// TODO: send 'nextWordCount' via a shared_ptr instead of pass it in and return it
nextWordCount mapWord (nextWordCount p, std::string key, std::string value) {
  std::cout << "Mapping KEY: " << key << " to VALUE: " << value << "\n";
  if (p.count(key)) {
    if (p[key].count(value)) {
      ++p[key][value];
    } else {
      p[key][value] = 1;
    }
  } else {
    p[key] = std::map<std::string, int> {{value, 1}};
  }
  return p;
}

nextWordProb calcFile (std::string filename) {
  std::ifstream fin (filename.c_str());
  if (!fin.is_open()) {
    std::cout << "File " << filename << " could not be opened\n";
    nextWordProb t;
    return t;
  }

  std::string cur, next;
  nextWordCount nwc;

  fin >> cur;
  while (fin >> next) {
    nwc = mapWord(nwc, cur, next);
    cur = next;
  }

  nextWordProb nwp = calcNextProb(nwc);

  return nwp;
}

// GENERATION NEW RANDOM TEXT
// given length and bias
// sudo:
//  cur = random key in bias
//  for size desired
//    next = possible_next_keys(cur)
//    write(cur)
//    cur = next;

std::string randomKey (nextWordProb nwp) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, nwp.size() -1);
  int jump = dis(gen);
  auto iter = nwp.begin();
  std::advance(iter, jump);
  return iter->first;
}

// REF: http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
std::string findNext (nextWordProb nwp, std::string cur) {
  if (!nwp.count(cur)) {
    std::cout << "ERROR: couldn't find " << cur << " in nwp\n";
    std::string t;
    return t;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(nwp[cur].second.begin(), nwp[cur].second.end());

  return nwp[cur].first[d(gen)];
}

nextWordProb calcNextProb (nextWordCount nwc) {
  nextWordProb nwp;

  for (auto nw : nwc) {
    std::vector<std::string> words;
    std::vector<int> weights;
    for (auto w : nw.second) {
      words.push_back(w.first);
      weights.push_back(w.second);
    }
    nwp[nw.first] = std::make_pair(words, weights);
  }
  return nwp;
}

void generateText (std::string filename, nextWordProb nwp, unsigned long long int size) {
  std::string cur = randomKey(nwp);
  std::ofstream fout (filename);
  fout << cur << " ";
  for (unsigned long long int i = 0; i < size; ++i) {
    cur = findNext(nwp, cur);
    fout << cur << " ";
  }
  fout.close();
}
