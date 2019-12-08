#include <vector>
#include "common.h"
#include <string>

std::vector<Node> loadFromFile(std::string filename, stepParams *params, int *num_iterations);

bool saveToFile(std::vector<Node> data, std::string filename);
