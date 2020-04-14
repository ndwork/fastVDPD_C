//
//  main.c
//  cimpl
//
//  Created by Nicholas Dwork starting on 9/18/16.
//  Copyright © 2016 Nicholas Dwork.
//

#include <stdlib.h>

#include "fastVDPD.h"



int main(int argc, const char * argv[]) {

  srand( 5 );  // Seed the random number generator

  //cimpl_runTests();

  runFastVDPD();

  return 0;
}
