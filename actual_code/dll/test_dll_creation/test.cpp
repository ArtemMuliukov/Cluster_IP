#define _SCL_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#define DllImport   __declspec( dllimport )  
#define DllExport   __declspec( dllexport ) 

#include <vector>
#include <list>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <iostream>


using namespace std;

extern int* tree_starts_saved = 0;
extern double* tree_saved = 0;
extern int m_size = 0;
extern int h_size = 0;
extern int tr_size = 0;

// struct contain the point data
struct point {
	vector<double> hiden; // Real object parameters
	vector<double> measured; // Mesured point parameters
	int sorted_param = 0; // advansed for build the tree
	friend bool operator<(point& a, point& b)
	{
		return a.measured[a.sorted_param] < b.measured[b.sorted_param];
	}
};

// vertex of hierarhical tree
struct vertex {
	unsigned int data_start = 0, data_end = 0;
	double R = 0;
	vector<double> hiden_center;
	vector<double> measured_center;
};

// transform point data to vector
vector<double> transform_point(double* point, int size) {
	vector<double> arr(size);
	for (int i = 0; i < size; ++i) {
		arr[i] = point[i];
	}
	return arr;
}

vector<point>* transform_mesured_data(double* points, int number_points,
	int hiden_size, int mesured_size) {
	vector<point>* res = new vector<point>(number_points);
	for (int i = 0; i < number_points; ++i) {
		(*res)[i].hiden = transform_point(points + i * (hiden_size + mesured_size), hiden_size);
		(*res)[i].measured = transform_point(points + hiden_size + i * (hiden_size + mesured_size), mesured_size);
	}
	return res;
}


// build false binary tree and save data in tree. Tree should be allocated before call [start; end)
int build(int start, int end, vector<point>* data, vector<vertex>* tree, int step) {

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
				* ((*tree)[step].measured_center[j] - (*data)[i].measured[j]);
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
	//if (end - start > 10000)
	//	printf("%d %d %d\n", start, end, dim);

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
		dist += (a[j] - b[j]) * (a[j] - b[j]);
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
//ratio_perm - the ration of R and sit dermited to reject, othercase we continue to split
//is_buid_tree - do we need the heap building
// return - num of probs (K)
DllExport int c_probs(int number_points, int hiden_size, int mesured_size,
	double* measure, double* centers, double* indicatrix, double* dist_s, int* weights, double* closest,
	double ratio_perm, int is_buid_tree, bool load, int max_clust_num) {

	//load if needed
	if (load) {
		FILE* fout = fopen("built_tree_data.bin", "rb");
		FILE* fconf = fopen("common_info.conf", "rt");
		long number_points;
		fscanf(fconf, "%ld", &number_points);
		fscanf(fconf, "%d", &hiden_size);
		fscanf(fconf, "%d", &mesured_size);

		m_size = mesured_size;
		h_size = hiden_size;
		tr_size = 2 << (int(logb(number_points)) + 1);
		tree_starts_saved = new int[tr_size * 2];
		tree_saved = new double[(1 + hiden_size + mesured_size) * tr_size];

		fread(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fread(tree_saved, sizeof(double), tr_size * (hiden_size + mesured_size + 1), fout);
		fclose(fconf);
		fclose(fout);
		return 0;
	}

	// asserts
	if (hiden_size != h_size)
		return -1;
	if (mesured_size != m_size)
		return -2;
	if (tree_saved == 0)
		return -3;

	//ofstream myfile;
	//myfile.open("debugdata.txt", std::ofstream::out | std::ofstream::app);

	vector<clust> cl_stack;// stack of clusters for checking
	clust first = { 0, count_dist(tree_saved + 1 + hiden_size, measure, mesured_size) };
	cl_stack.push_back(first);
	if (is_buid_tree)
		make_heap(cl_stack.begin(), cl_stack.end());


	int probs_counter = 0;
	double upper_bound_dist = DBL_MAX;
	double smallest_dist = DBL_MAX;
	//vector<int> nums_of_p_in_clust;
	int comps_num = 0;
	while (cl_stack.size() > 0) {
		comps_num++;
		int split = 1;
		clust cur;
		if (is_buid_tree)
			cur = cl_stack[0];
		else
			cur = cl_stack[cl_stack.size() - 1];
		double dist = count_dist(tree_saved + (hiden_size + mesured_size + 1) * cur.num + 1 + hiden_size,
			measure, mesured_size);

		// Work with this element
		//******

		// do it if cluster is to be thrown
		if ((upper_bound_dist < dist - (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0] &&
			ratio_perm*dist > (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0]) ||
			tree_starts_saved[2 * (cur.num) + 1] - tree_starts_saved[2 * (cur.num)] == 1) { // here check of splitting
			if (max_clust_num <= probs_counter)
				return -4;
			split = 0;
			dist_s[probs_counter] = dist*dist;
			weights[probs_counter] = tree_starts_saved[2 * cur.num + 1] - tree_starts_saved[2 * cur.num];
			//currently i pushed dist, not the probs
			copy(tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1,
				tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1 + hiden_size,
				centers + probs_counter * hiden_size);
			copy(tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1 + hiden_size,
				tree_saved + (cur.num + 1) * (hiden_size + mesured_size + 1),
				indicatrix + probs_counter * hiden_size);
			probs_counter++;
			/*myfile << probs_counter << ' ' << cur.num << ' ' << cur.dist << " " <<
				upper_bound_dist << " " << cl_stack.size() << " " << comps_num << std::endl;*/
				//nums_of_p_in_clust.push_back(tree_starts_saved[2 * cur.num + 1] - tree_starts_saved[2 * cur.num]);
				// smallest dist chandge
		}

		// new upped bound
		if (upper_bound_dist > dist + (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0])
			upper_bound_dist = dist + (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0];
		//*****
		if (smallest_dist > dist + *(tree_saved + cur.num * (hiden_size + mesured_size + 1))) {
			smallest_dist = dist + *(tree_saved + cur.num * (hiden_size + mesured_size + 1));
			copy(tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1,
				tree_saved + (cur.num + 1) * (hiden_size + mesured_size + 1),
				closest);
		}
		//delete it from stack
		if (is_buid_tree)
			pop_heap(cl_stack.begin(), cl_stack.end());
		cl_stack.pop_back();
		//put new elements
		if (split) {
			clust sec = { 2 * cur.num + 1, count_dist(tree_saved + (2 * cur.num + 1) * (1 + hiden_size + mesured_size)
				+ 1 + hiden_size, measure, mesured_size) /*- tree[(2 * cur.num + 1)*(1 + hiden_size + mesured_size)]*/ };
			clust third = { 2 * cur.num + 2, count_dist(tree_saved + (2 * cur.num + 2) * (1 + hiden_size + mesured_size)
				+ 1 + hiden_size, measure, mesured_size) /*- tree[(2 * cur.num + 2)*(1 + hiden_size + mesured_size)]*/ };
			if (is_buid_tree) {
				cl_stack.push_back(sec); push_heap(cl_stack.begin(), cl_stack.end());
				cl_stack.push_back(third); push_heap(cl_stack.begin(), cl_stack.end());
			}
			else {
				if (sec.dist <= third.dist) {
					cl_stack.push_back(third);
					cl_stack.push_back(sec);
				}
				if (sec.dist > third.dist) {
					cl_stack.push_back(sec);
					cl_stack.push_back(third);
				}
			}
		}
	}

	return probs_counter;
	//printf(" %d ", comps_num);
	//return comps_num;
}


DllExport int c_probs_fast(int number_points, int hiden_size, int mesured_size,
	double* measure, double* centers, double* dist_s, int* weights, double* closest,
	double ratio_perm, int is_buid_tree, bool load, int max_clust_num) {

	//load if needed
	if (load) {
		FILE* fout = fopen("built_tree_data.bin", "rb");
		FILE* fconf = fopen("common_info.conf", "rt");
		long number_points;
		fscanf(fconf, "%ld", &number_points);
		fscanf(fconf, "%d", &hiden_size);
		fscanf(fconf, "%d", &mesured_size);

		m_size = mesured_size;
		h_size = hiden_size;
		tr_size = 2 << (int(logb(number_points)) + 1);
		tree_starts_saved = new int[tr_size * 2];
		tree_saved = new double[(1 + hiden_size + mesured_size) * tr_size];

		fread(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fread(tree_saved, sizeof(double), tr_size * (hiden_size + mesured_size + 1), fout);
		fclose(fconf);
		fclose(fout);
		return 0;
	}

	// asserts
	if (hiden_size != h_size)
		return -1;
	if (mesured_size != m_size)
		return -2;
	if (tree_saved == 0)
		return -3;

	//ofstream myfile;
	//myfile.open("debugdata.txt", std::ofstream::out | std::ofstream::app);

	vector<clust> cl_stack;// stack of clusters for checking
	clust first = { 0, count_dist(tree_saved + 1 + hiden_size, measure, mesured_size) };
	cl_stack.push_back(first);
	if (is_buid_tree)
		make_heap(cl_stack.begin(), cl_stack.end());

	//list<int> cl_stack;// stack of clusters for checking, old variant
	//cl_stack.push_back(0);
	int probs_counter = 0;
	double upper_bound_dist = DBL_MAX;
	double smallest_dist = DBL_MAX;
	//vector<int> nums_of_p_in_clust;
	int comps_num = 0;
	while (cl_stack.size() > 0) {
		comps_num++;
		int split = 1;
		clust cur;
		if (is_buid_tree)
			cur = cl_stack[0];
		else
			cur = cl_stack[cl_stack.size() - 1];
		double dist = count_dist(tree_saved + (hiden_size + mesured_size + 1) * cur.num + 1 + hiden_size,
			measure, mesured_size);

		// Work with this element
		//******

		// do it if cluster is to be thrown
		if ((upper_bound_dist < dist - (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0] &&
			ratio_perm * dist >(tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0]) ||
			tree_starts_saved[2 * (cur.num) + 1] - tree_starts_saved[2 * (cur.num)] == 1) { // here check of splitting
			if (max_clust_num <= probs_counter)
				return -4;
			split = 0;
			dist_s[probs_counter] = dist*dist;
			weights[probs_counter] = tree_starts_saved[2 * cur.num + 1] - tree_starts_saved[2 * cur.num];
			//currently i pushed dist, not the probs
			copy(tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1,
				tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1 + hiden_size,
				centers + probs_counter * hiden_size);
			probs_counter++;
			/*myfile << probs_counter << ' ' << cur.num << ' ' << cur.dist << " " <<
				upper_bound_dist << " " << cl_stack.size() << " " << comps_num << std::endl;*/
				//nums_of_p_in_clust.push_back(tree_starts_saved[2 * cur.num + 1] - tree_starts_saved[2 * cur.num]);
				// smallest dist chandge
		}

		// new upped bound
		if (upper_bound_dist > dist + (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0])
			upper_bound_dist = dist + (tree_saved + (hiden_size + mesured_size + 1) * cur.num)[0];
		//*****
		if (smallest_dist > dist + *(tree_saved + cur.num * (hiden_size + mesured_size + 1))) {
			smallest_dist = dist + *(tree_saved + cur.num * (hiden_size + mesured_size + 1));
			copy(tree_saved + cur.num * (hiden_size + mesured_size + 1) + 1,
				tree_saved + (cur.num + 1) * (hiden_size + mesured_size + 1),
				closest);
		}
		//delete it from stack
		if (is_buid_tree)
			pop_heap(cl_stack.begin(), cl_stack.end());
		cl_stack.pop_back();
		//put new elements
		if (split) {
			clust sec = { 2 * cur.num + 1, count_dist(tree_saved + (2 * cur.num + 1) * (1 + hiden_size + mesured_size)
				+ 1 + hiden_size, measure, mesured_size) /*- tree[(2 * cur.num + 1)*(1 + hiden_size + mesured_size)]*/ };
			clust third = { 2 * cur.num + 2, count_dist(tree_saved + (2 * cur.num + 2) * (1 + hiden_size + mesured_size)
				+ 1 + hiden_size, measure, mesured_size) /*- tree[(2 * cur.num + 2)*(1 + hiden_size + mesured_size)]*/ };
			if (is_buid_tree) {
				cl_stack.push_back(sec); push_heap(cl_stack.begin(), cl_stack.end());
				cl_stack.push_back(third); push_heap(cl_stack.begin(), cl_stack.end());
			}
			else {
				if (sec.dist <= third.dist) {
					cl_stack.push_back(third);
					cl_stack.push_back(sec);
				}
				if (sec.dist > third.dist) {
					cl_stack.push_back(sec);
					cl_stack.push_back(third);
				}
			}
		}
	}

	return probs_counter;
	//printf(" %d ", comps_num);
	//return comps_num;
}

// function work with allocated arrays
// data input data in format {hiden..., musured...}^N
// tree_starts - indexes in data {start, end}^M, M - number of knots (as max, twice bigger than N)
// tree - points in tree. fromat: {R, hiden..., musured...}^M 
DllExport int32_t build_tree(double* data, int32_t number_points,
	int32_t hiden_size, int32_t mesured_size, int32_t save, int32_t use_saved) {

	if (use_saved) {
		FILE* fout = fopen("built_tree_data.bin", "rb");
		FILE* fconf = fopen("common_info.conf", "rt");
		fscanf(fconf, "%ld", &number_points);
		fscanf(fconf, "%d", &hiden_size);
		fscanf(fconf, "%d", &mesured_size);

		m_size = mesured_size;
		h_size = hiden_size;
		tr_size = 2 << (int(logb(number_points)) + 1);
		tree_starts_saved = new int[tr_size * 2];
		tree_saved = new double[(1 + hiden_size + mesured_size) * tr_size];

		fread(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fread(tree_saved, sizeof(double), tr_size * (hiden_size + mesured_size + 1), fout);
		fclose(fconf);
		fclose(fout);
		return 1;
	}

	tr_size = 2 << (int(logb(number_points)) + 1);

	vector<point>* data_vect = transform_mesured_data(data, number_points, hiden_size, mesured_size);
	vector<vertex>* tree_vect = new vector<vertex>(tr_size);

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

	tree_starts_saved = new int[tr_size * 2];
	tree_saved = new double[(1 + hiden_size + mesured_size) * tr_size];

	m_size = mesured_size;
	h_size = hiden_size;

	for (int i = 0; i < tree_vect->size(); ++i) {
		tree_starts_saved[2 * i] = (*tree_vect)[i].data_start;
		tree_starts_saved[2 * i + 1] = (*tree_vect)[i].data_end;

		tree_saved[i * (1 + hiden_size + mesured_size)] = (*tree_vect)[i].R;
		for (int j = 0; j < (*tree_vect)[i].hiden_center.size(); ++j) {
			tree_saved[i * (1 + hiden_size + mesured_size) + 1 + j] = (*tree_vect)[i].hiden_center[j];
		}
		for (int j = 0; j < (*tree_vect)[i].measured_center.size(); ++j) {
			tree_saved[i * (1 + hiden_size + mesured_size) + 1 + hiden_size + j] = (*tree_vect)[i].measured_center[j];
		}
	}

	if (save) {
		FILE* fout = fopen("built_tree_data.bin", "wb");
		FILE* fconf = fopen("common_info.conf", "wt");
		fprintf(fconf, "%ld ", number_points);
		fprintf(fconf, "%d ", hiden_size);
		fprintf(fconf, "%d ", mesured_size);
		//fwrite(data, sizeof(double), database_size * (hid_size + mes_size), fout);
		fwrite(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fwrite(tree_saved, sizeof(double), tr_size * (hiden_size + mesured_size + 1), fout);
		fclose(fconf);
		fclose(fout);	
	}
	return 1;
}

double maltsev_norm_func(double num) {
	return exp(-2 * log(num / 54) * log(num / 54)) / num;
}


//////////////////////////////////////////////////
//lower the tests
////////////////////////////////////////////////////

// Fuction for test of saving tree data
DllExport int get_msize(int var) {
	return m_size;
}


double* load_binary_data(const char* name, int* num, int* hiden, int* measured) {
	FILE* f = fopen(name, "rb");
	fread(num, sizeof(int), 1, f);
	fread(hiden, sizeof(int), 1, f);
	fread(measured, sizeof(int), 1, f);
	int double_nums = *num * (*hiden + *measured);
	double* data = new double[double_nums];
	auto r = fread(data, sizeof(double), double_nums, f);
	return data;
}

double* load_data(int string_size, const char* name, int* num_of_elements, int norm) {
	double* res = new double[string_size * (*num_of_elements)];
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
	double* data = load_binary_data("databasebytes_new", &database_size, &hid_size, &mes_size);

	int upperbound_tree = 2 << (int(logb(database_size)) + 1);

	build_tree(data, database_size, hid_size, mes_size, true, 0);

	return;
}
