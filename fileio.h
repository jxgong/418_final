#ifndef FILEIO_H
#define FILEIO_H

#include <vector>
#include "common.h"
#include <string>

std::vector<Node> loadFromFile(std::string filename);

bool saveToFile(std::vector<Node> data, std::string filename);
#endif
