#include <vector>
#include "common.h"
#include <string>

#ifndef FILEIO_H
#define FILEIO_H

std::vector<Node> loadFromFile(std::string filename, stepParams *params, int *num_iterations);

bool saveToFile(std::vector<Node> data, std::string filename);
#endif
