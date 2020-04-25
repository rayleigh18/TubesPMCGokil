#include <stdio.h>
#include <stdlib.h>
#include "configuration.h"

int main() {
  component *com;
  int len;
  char typ;
  double val;
  int n1;
  int n2;
  double t;
//  *com = (component**)malloc(sizeof((component*)*5);
 // com = ;
//  **com.node1=0;
 // **com.node2=0;
 // **com.value=0;
//  **com.type=' ';

  typ = 'L';
  val = 5.5;
  n1 = 2;
  n2 = 3;
  t = 3;
initializeComponentArray(&com);
addComponent(&com, &len, typ, val, n1, n2, t);
printComponents(com, len, t);
destroyComponentArray(&com);
}
