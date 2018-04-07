#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace std;

double maltsev_norm_func(double num) {
	return exp(-2 * log(num / 54)* log(num / 54)) / num;
}

int transform_file() {
	FILE* inf;
	FILE* outf;
	inf = fopen("new_db_61.data", "rt");
	outf = fopen("databasebytes_new", "wb");
	if (inf)
		cout << "in is opened\n";
	if (outf)
		cout << "out is opened\n";
	int counter = 0 ;
	int bound = 500000 * 65;

	int num = 198877;//493739;
	fwrite(&num, sizeof(int), 1, outf); //write the number of elements
	num = 4;
	fwrite(&num, sizeof(int), 1, outf); //write the number of hides points
	num = 61;
	fwrite(&num, sizeof(int), 1, outf); //write the number of measured points

	double x;
	while (fscanf(inf, "%lf", &x) == 1) {
		//if (counter % 65 >=4)
			//x *= maltsev_norm_func(counter % 65 + 6);//normation the data, 6 = 10(grad) - 4(num of hiden)
		fwrite(&x, sizeof(double), 1, outf);//пишем по 8 байтов double в файл
		counter++;
		if (counter % 325000 == 0)
			printf("%d%% \n", counter / 325000);
	}
	printf("Read %d elements", counter);
	fclose(outf);
	return 9;
}

void check_reading() {
	FILE * infile = fopen("databasebytes_new", "rb");
	int num = 199877;//493739;
	double * data = new double[num * 65];
	auto r = fread(data, sizeof(double), num * 65, infile);
	printf("olala");
	5 + 6;
}

// этот код поможет мне перевести файл из формата чисел в текстовом файле в формат файла из байтов-чисел
int main() {
	transform_file();
	check_reading();
}