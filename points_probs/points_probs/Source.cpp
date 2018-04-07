#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS

#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <time.h>


using namespace std;

// struct contain the point data
struct point {
	vector<double> hiden; // Real object parameters
	vector<double> measured; // Mesured point parameters
	int sorted_param; // advansed for build the tree
	friend bool operator<(point& a, point& b)
	{
		return a.measured[a.sorted_param] < b.measured[b.sorted_param];
	}
};

// vertex of hierarhical tree
struct vertex {
	unsigned int data_start, data_end;
	double R;
	vector<double> hiden_center;
	vector<double> measured_center;
};

// transform point data to vector
vector<double> transform_point(double *point, int size) {
	vector<double> arr(size);
	for (int i = 0; i < size; ++i) {
		arr[i] = point[i];
	}
	return arr;
}

vector<point>* transform_mesured_data(double *points, int number_points,
	int hiden_size, int mesured_size) {
	vector<point> *res = new vector<point>(number_points);
	for (int i = 0; i < number_points; ++i) {
		(*res)[i].hiden = transform_point(points + i * (hiden_size + mesured_size), hiden_size);
		(*res)[i].measured = transform_point(points + hiden_size + i * (hiden_size + mesured_size), mesured_size);
	}
	return res;
}


// build false binary tree and save data in tree. Tree should be allocated before call [start; end)
int build(int start, int end, vector<point> *data, vector<vertex> *tree, int step) {

	//make vertex
	(*tree)[step].data_start = start;
	(*tree)[step].data_end = end;
	for (int i = 0; i < (*data)[0].hiden.size(); ++i) {
		(*tree)[step].hiden_center.push_back(0); // new coordinate
		for (int j = start; j < end; ++j) {
			(*tree)[step].hiden_center[i] += (*data)[j].hiden[i]; //sum around start-end
		}
		(*tree)[step].hiden_center[i] /= end - start;//normalize
	}
	for (int i = 0; i < (*data)[0].measured.size(); ++i) {
		(*tree)[step].measured_center.push_back(0);
		for (int j = start; j < end; ++j) {
			(*tree)[step].measured_center[i] += (*data)[j].measured[i];
		}
		(*tree)[step].measured_center[i] /= end - start;
	}
	double max_rad = 0;
	for (int i = start; i < end; ++i) {
		double rad = 0;
		for (int j = 0; j < (*tree)[step].measured_center.size(); ++j) {
			rad += ((*tree)[step].measured_center[j] - (*data)[i].measured[j])
				*((*tree)[step].measured_center[j] - (*data)[i].measured[j]);
		}
		if (rad > max_rad)
			max_rad = rad;
	}
	(*tree)[step].R = sqrt(max_rad);
	//////////////////

	// searching biggest dimention-span {dim}
	double max_gap = 0;
	int dim = 0;
	for (int i = 0; i < (*data)[start].measured.size(); ++i) {
		double min = (*data)[start].measured[i];
		double max = (*data)[start].measured[i];
		for (int j = start; j < end; ++j) {
			if ((*data)[j].measured[i] < min)
				min = (*data)[j].measured[i];
			if ((*data)[j].measured[i] > max)
				max = (*data)[j].measured[i];
		}
		if (max_gap < max - min) {
			max_gap = max - min;
			dim = i;
		}
	}
	//debag
	if (end - start > 10000)
		printf("%d %d %d\n", start, end, dim);

	for (int i = start; i < end; ++i) {
		(*data)[i].sorted_param = dim;
	}

	// medianth-stat partical sort and recursive call
	int mid = (end + start) / 2;
	if (end - start >= 2) {
		nth_element((*data).begin() + start, (*data).begin() + mid, (*data).begin() + end);
		build(start, mid, data, tree, 2 * step + 1);
		build(mid, end, data, tree, 2 * step + 2);
	}
	return 0;
}

struct clust {
	int num;
	float dist;
	friend bool operator<(clust a, clust b)
	{
		return a.dist > b.dist;
	}
};

double count_dist(double* a, double* b, int num) {
	double dist = 0;
	for (int j = 0; j < num; ++j) {
		dist += (a[j] - b[j])*(a[j] - b[j]);
	}
	dist = sqrt(dist);
	return dist;
}


//////////////////////////////////////////////////////////////////////////////
// next libs work in dll
////////////////////////////////////////////////////////////



// function work with allocated arrays
// data input data in format {hiden..., musured...}^N
// tree_starts - indexes in data {start, end}^M, M - number of knots (as max, twice bigger than N)
// tree - points in tree. fromat: {R, hiden..., musured...}^M 
// measure - measured points in fromat {mearured...}
// above array are ready, next is result of work
// centers - centers of probs {hiden...}^K, K - number of return probs. less than num of points
// probs - {1}^K
// return - num of probs (K)
int c_probs(int *tree_starts, double *tree, int number_points,
	int hiden_size, int mesured_size, double * measure,	double * centers, double * probs) {


	vector<clust> cl_stack;// stack of clusters for checking
	clust first = { 0,count_dist(tree + 1 + hiden_size, measure, mesured_size) - tree[0] };
	cl_stack.push_back(first);
	make_heap(cl_stack.begin(), cl_stack.end());

	//list<int> cl_stack;// stack of clusters for checking, old variant
	//cl_stack.push_back(0);
	int probs_counter = 0;
	double upper_bound_dist = DBL_MAX;
	vector<int> nums_of_p_in_clust;
	int comps_num = 0;
	while (cl_stack.size() > 0) {
		comps_num++;
		int split = 1;
		clust cur = cl_stack[0];
		double dist = count_dist(tree + (hiden_size + mesured_size + 1)*cur.num + 1 + hiden_size, measure, mesured_size);
		// Work with this element
		//******

		// do it if cluster is throw
		if (upper_bound_dist < dist - (tree + (hiden_size + mesured_size + 1)*cur.num)[0] ||
			tree_starts[2 * (cur.num) + 1] - tree_starts[2 * (cur.num)] == 1) { // here check of splitting

			split = 0;
			probs[probs_counter] = dist; 
			//currently i pushed dist but it must be probs!!!
			copy(tree + cur.num*(hiden_size + mesured_size + 1) + 1,
				tree + (cur.num+1) * (hiden_size + mesured_size + 1) + 1,
				centers + probs_counter * (hiden_size + mesured_size));
			probs_counter++;
			nums_of_p_in_clust.push_back(tree_starts[2 * cur.num + 1] - tree_starts[2 * cur.num]);
		}

		// new upped bound
		if (upper_bound_dist > dist + (tree + (hiden_size + mesured_size + 1)*cur.num)[0])
			upper_bound_dist = dist + (tree + (hiden_size + mesured_size + 1)*cur.num)[0];
		//*****
		//delete it from stack
		pop_heap(cl_stack.begin(), cl_stack.end()); cl_stack.pop_back();
		//put new elements
		if (split) {
			clust sec = { 2 * cur.num + 1,count_dist(tree + (2 * cur.num + 1)*(1 + hiden_size + mesured_size)
				+ 1 + hiden_size, measure, mesured_size) - tree[(2 * cur.num + 1)*(1 + hiden_size + mesured_size)] };
			clust third = { 2 * cur.num + 2,count_dist(tree + (2 * cur.num + 2)*(1 + hiden_size + mesured_size)
				+ 1 + hiden_size, measure, mesured_size) - tree[(2 * cur.num + 2)*(1 + hiden_size + mesured_size)] };
			cl_stack.push_back(sec); push_heap(cl_stack.begin(), cl_stack.end());
			cl_stack.push_back(third); push_heap(cl_stack.begin(), cl_stack.end());
		}
	}
	//double sum_pr = 0;
	//double k_eff = 20; // here should be func of calc k_eff
	////summarize probs
	//for (int i = 0; i < probs_counter; i++) {
	//	probs[i] = pow(probs[i], -k_eff / 2)*nums_of_p_in_clust[i];
	//	sum_pr += probs[i];
	//}
	//for (int i = 0; i < probs_counter; i++)
	//	probs[i] /= sum_pr;

	//printf(" %d ", comps_num);
	return probs_counter;
	//return comps_num;
}


// function work with allocated arrays
// data input data in format {hiden..., musured...}^N
// tree_starts - indexes in data {start, end}^M, M - number of knots (as max, twice bigger than N)
// tree - points in tree. fromat: {R, hiden..., musured...}^M 
int build_tree(double *data, int *tree_starts, double *tree, int number_points,
	int hiden_size, int mesured_size) {
	vector<point> *data_vect = transform_mesured_data(data, number_points, hiden_size, mesured_size);
	vector<vertex> *tree_vect = new vector<vertex>(2 << (int(logb(number_points)) + 1));

	build(0, number_points, data_vect, tree_vect, 0);

	// saving database/ not used
	/*for (int i = 0; i < data_vect->size(); ++i) {
		memcpy(&(data[i * (hiden_size + mesured_size)]),
			(*data_vect)[i].hiden.data(),
			hiden_size * sizeof(double));
		memcpy(&(data[i * (hiden_size + mesured_size) + hiden_size]),
			(*data_vect)[i].measured.data(),
			mesured_size * sizeof(double));
	}*/
	for (int i = 0; i < tree_vect->size(); ++i) {
		tree_starts[2 * i] = (*tree_vect)[i].data_start;
		tree_starts[2 * i + 1] = (*tree_vect)[i].data_end;

		tree[i * (1 + hiden_size + mesured_size)] = (*tree_vect)[i].R;
		for (int j = 0; j < (*tree_vect)[i].hiden_center.size(); ++j) {
			tree[i * (1 + hiden_size + mesured_size) + 1 + j] = (*tree_vect)[i].hiden_center[j];
		}
		for (int j = 0; j < (*tree_vect)[i].measured_center.size(); ++j) {
			tree[i * (1 + hiden_size + mesured_size) + 1 + hiden_size + j] = (*tree_vect)[i].measured_center[j];
		}
	}
	return 0;
}

double maltsev_norm_func(double num) {
	return exp(-2 * log(num / 54)* log(num / 54)) / num;
}

//////////////////////////////////////////////////
//lower the tests
////////////////////////////////////////////////////


double *load_binary_data(const char* name, int* num, int *hiden, int *measured) {
	FILE *f = fopen(name, "rb");
	fread(num, sizeof(int), 1, f);
	fread(hiden, sizeof(int), 1, f);
	fread(measured, sizeof(int), 1, f);
	int double_nums = *num * (*hiden + *measured);
	double * data = new double[double_nums];
	auto r = fread(data, sizeof(double), double_nums, f);
	return data;
}

double *load_data(int string_size, const char* name, int *num_of_elements , int norm) {
	double *res = new double[string_size * (*num_of_elements)];
	FILE* file = fopen(name, "rt");
	int counter = 0;
	double num;
	while (fscanf(file, "%lf", res + counter) == 1) {
		if (norm)
			res[counter] *= maltsev_norm_func(counter % 65 - 4 + 10);
		counter++;
	}
	*num_of_elements = counter / string_size;
	fclose(file);
	return res;
}


void build_and_save_tree() {
	int hid_size = 4;
	int mes_size = 61;
	int database_size = 300000;

	// read full database, normed not in ahother programm
	double *data = load_binary_data("databasebytes_new", &database_size, &hid_size, &mes_size);

	int upperbound_tree = 2 << (int(logb(database_size)) + 1);
	//build testing tree
	int * tree_starts = new int[upperbound_tree * 2];
	double * tree = new double[upperbound_tree *(hid_size + mes_size + 1)];

	build_tree(data, tree_starts, tree, database_size, hid_size, mes_size);

	//save tree build
	FILE * fout = fopen("built_tree_data", "wb");
	FILE * fconf = fopen("common_info.conf", "wt");
	fprintf(fconf, "%d ", database_size);
	fprintf(fconf, "%d ", hid_size);
	fprintf(fconf, "%d ", mes_size);
	//fwrite(data, sizeof(double), database_size * (hid_size + mes_size), fout);
	fwrite(tree_starts, sizeof(double), upperbound_tree, fout);
	fwrite(tree, sizeof(double), upperbound_tree * (hid_size + mes_size + 1), fout);
	fclose(fconf);
	fclose(fout);

	return;
}

void test_f_find_closest() {

	int hid_size = 4;
	int mes_size = 61;
	int database_size = 300000;

	FILE * fconf = fopen("common_info.conf", "rt");
	fscanf(fconf, "%d", &database_size);
	fscanf(fconf, "%d", &hid_size);
	fscanf(fconf, "%d", &mes_size);

	int upperbound_tree = 2 << (int(logb(database_size)) + 1);

	//double *data = new double[database_size * (hid_size + mes_size)];
	int * tree_starts = new int[upperbound_tree * 2];
	double * tree = new double[upperbound_tree *(hid_size + mes_size + 1)];

	FILE * fin = fopen("built_tree_data", "rb");
	//fread(data, sizeof(double), database_size * (hid_size + mes_size), fin);
	fread(tree_starts, sizeof(double), upperbound_tree, fin);
	fread(tree, sizeof(double), upperbound_tree * (hid_size + mes_size + 1), fin);
	fclose(fin);
	fclose(fconf);

	int test_size = 1000;
	double *data_mod = load_data(hid_size + mes_size, "database_txt/test_1000.txt",&test_size, 0); // read model database
	double *data_exp = load_data(mes_size, "database_txt/experiment_1000.txt", &test_size, 0); // read exp database

	auto t = clock();
	// create file for outs
	//FILE * file = fopen("results.txt", "wt");
	//make place for answers	
	double *centers = new double[database_size * (hid_size + mes_size)];
	double *probs = new double[database_size];
	int sum = 0;
	// calc probs for each indicatrix
	for (int i = 0; i < test_size; i++) {
		int num_probs = c_probs(tree_starts, tree, database_size, hid_size, mes_size,
			data_exp + i * mes_size, centers, probs);
		//printf("%f %%  %d\n", float(i) / 10, num_probs);
		sum += num_probs;
		/*double max = 100;
		int ind_max = -1;
		for (int k = 0; k < num_probs; k++)
			if (max > probs[k]) {
				max = probs[k];
				ind_max = k;
			}
		fprintf(file, "%f ", max);
		for (int k = 0; k < hid_size + mes_size; ++k) {
			fprintf(file, "%f ", centers[ind_max*(hid_size + mes_size) + k]);
		}
		fprintf(file, "\n");*/
	}
	printf("\n Time gone: %f sec.\n", (clock() - t) / 1000.0);
	printf(" %d\n", sum);
	return;
}

int main() {
	printf("Start");
	//build_and_save_tree();
	test_f_find_closest();
	printf("Finish?");
	getchar();
	return 0;
}