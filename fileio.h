#include <vector>
#include "common.h"
#include <string>
#include <unordered_set>

#ifndef FILEIO_H
#define FILEIO_H

std::vector<Node> loadFromFile(std::string filename, stepParams *params, int *num_iterations, std::unordered_set<int> *sparks, int *sparkStart, int *sparkEnd);

bool saveToFile(std::vector<Node> data, std::string filename);
#endif
