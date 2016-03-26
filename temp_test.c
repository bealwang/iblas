#include <stdio.h>
#include <stdlib.h>
#include "mic/icpc.h"

// #define LEN 5

// __attribute__((target(mic))) void funcheck(int i){
// #ifdef __MIC__
// 	printf("INDEX on MIC:%d\n",i);
// #else 
// 	printf("INDEX on CPU:%d\n",i);
// #endif
// }

// __attribute__((target(mic))) void check(){
// 	float x = 2;
// 	float arr[LEN];
// 	int i;

// #pragma offload target(mic)
// 	for(i=0;i<LEN;++i){
// 		//arr[i] = i*3.0f/2.0f;
// 		funcheck(i);
// 	}
// }

int main(){
	check();
  	return 0;
}


