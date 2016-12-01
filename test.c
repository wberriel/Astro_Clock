#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{

	double jme = 2452930.312847;
	double A = 3341656.0;
	double B = 4.6692568;
	double C = 6283.07585;

	double result = 0;
	double cosResult = 0;
	double testa;
	double testb;
	//double testc;
	
	result = jme * C;
	cosResult = cos(result + B);
	testa = A * cosResult;
	testb = A * cos(B + C*jme);

	printf("result = %f cosResult = %f\n", result, cosResult);
	printf("TestA = %f, TestB = %f\n", testa, testb);
}

