#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <sys/stat.h> 
#include <sys/types.h> 


#include "config.h"
#include "function.h"
#include "parallel.h"

int main(int argc, char** argv) {
    Config conf(argc, argv);

    Function3D u(conf.getLx(), conf.getLy(), conf.getLz());

    Parallel q(conf, u);
    q.process();
    
    return 0; 
}