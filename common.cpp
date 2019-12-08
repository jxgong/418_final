#include "common.h"

StartupOptions parseOptions(int argc, char *argv[]){
    StartupOptions res;
    res.numIterations = 2;
    res.length = 10;
    res.width = 10;
    res.height = 1;
    res.inputFile = "input.txt";
    res.checkCorrectness = true;
    return res;
}
