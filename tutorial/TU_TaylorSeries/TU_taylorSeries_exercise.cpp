/* comment section

*/


#include "../../include/myNP_22100252.h"


int main(int argc, char* argv[])
{

	double x = PI / 3;
	//double x = 60;

	double S_N = 0;

	/*===== Select the function to call =====*/
	//S_N = sinTaylor(x);
	//S_N = sindTaylor(x);
	S_N = cosTaylor(x);
	//S_N = expTaylor(x);

	printf("\n\n");
	printf("=======================================\n");
	printf("    sin( %f[rad] ) Calculation   \n", x);
	printf("=======================================\n");
	printf("   -  My     result = %3.12f    \n", S_N);
	printf("   -  Math.h result = %3.12f    \n", sin(x));
	printf("   -  absolute err. = %3.12f    \n", S_N - sin(x));
	printf("=======================================\n");

	system("pause");
	return 0;
}


