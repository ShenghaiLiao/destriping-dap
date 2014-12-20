#include <stdio.h>


// parse the interval of integers from n0 to n1 with increment inc
// and find the number whose largest prime factor is the smallest
// of all the numbers parsed
int get_best_composite_number(int n0, int n1, int inc){
  int i, j, i1, j1, best_num = n0;
  int largest_prime = n0;
  
  for(i=n0; i<n1; i+=inc){
    // factorize i
    i1 = i; 
    for(j=2; j*j<i; j++){
      while( (i1%j)==0 ){ i1 /= j; }
      if(i1==1) break;
    }
    if(i1>1) j1 = i1; else j1 = j;
    // now j1 is the largest prime factor of i
    // update largest prime and "best number" if j1 is smaller
    if(largest_prime>j1) { largest_prime = j1; best_num = i; }
  }
  return best_num;
}



int main(int argc, char** argv) {

  int n = 1000;
  int inc = 2;
  float f1 = 1.01;
  float f2 = 1.1;

  if(argc>1) n  = atoi(argv[1]);
  if(argc>2) f1 = atof(argv[2]);

  int n1 = (n*f1);
  
  int i1, j, nn, power; 
 
  for(n=500; n<10000; n+=2) {
    nn = get_best_composite_number( n + 16, (int) (n*1.1), 2 );
    printf("You say %i   I say %i  = 1", n, nn);
    i1 = nn;
    for(j=2; j*j<nn; j++){
      power = 0;
      while( (i1%j)==0 ){ i1 /= j; power++; /* printf(" * %i", j); */ }
      if(power>0) { printf(" * %i^%i", j, power); }
      if(i1==1) break;
    }
    printf("\n");
    
  };

  return 0;
}
