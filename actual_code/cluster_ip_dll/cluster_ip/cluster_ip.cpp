/*
* CLUSTER_IP 0.3, under GPL-3.0 License.
*
* The solution is designed as C++ code in Microsoft Visual Studio. 

* The library is close to classical realisation of kd-tree, 
* with addition of flexible search, returning varying number of considered clusters 
* to enrich the solution (nearest element) with it's statistical estimates 
* (like mathematical expectation and standard deviation). 
* Also you could vary the depth of search to set appropriate compromise between speed and accuracy 
* for the statistical estimates (see the choice of ratio parameter in c_probs function).

* The code is compiled into a single dll file that has 2 key functions.

* "build_tree" - allows you to load an array of data into the dll, 
* which will be converted into a binary pseudo-tree and stored in the 
* internal memory of the DLL for subsequent operation of the program.

* "c_probs" - a function that finds the nearest element in the tree,
* as well as outputs the considered and discarded clusters (and distances to them) along the way, 
* which can be used for subsequent calculations (an example can be found in the folder with tests in LabVIew)

* The effectiveness of this implementation was tested on a real biological problem
* and presented at various conferences(e.g., ELS-XVII) and as a qualification master's thesis of Muliukov Artem.

*/

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
#include <string>

using namespace std;

// Variables used by both functions, they hold information about the constructed tree
extern int* tree_starts_saved = 0;
extern double* tree_saved = 0;
extern int m_size = 0;
extern int h_size = 0;
extern int tr_size = 0;

// The struct contains one point of tree for an intermediate processing of build_tree fuction
struct point {
	vector<double> hidden; // Real object parameters
	vector<double> measured; // measured point parameters (the measured signal)
	int sorted_param = 0; // the parameter number, along which the latest sort was performed (when building the tree)
	friend bool operator<(point& a, point& b)
	{
		return a.measured[a.sorted_param] < b.measured[b.sorted_param];
	}
};

// vertex of hierarhical tree
struct vertex {
	unsigned int data_start = 0, data_end = 0;
	double R = 0;
	vector<double> hidden_center;
	vector<double> measured_center;
};

// transform point data to vector (from array to vector)
vector<double> transform_point(double* point, int size) {
	vector<double> arr(size);
	for (int i = 0; i < size; ++i) {
		arr[i] = point[i];
	}
	return arr;
}

//array to object point
vector<point>* transform_measured_data(double* points, int number_points, int hidden_size, int measured_size) {
	vector<point>* res = new vector<point>(number_points);
	for (int i = 0; i < number_points; ++i) {
		(*res)[i].hidden = transform_point(points + i * (hidden_size + measured_size), hidden_size);
		(*res)[i].measured = transform_point(points + hidden_size + i * (hidden_size + measured_size), measured_size);
	}
	return res;
}

// build pseudo-binary tree and save data in tree. The tree memory MUST be allocated before call [start; end)
int build(int start, int end, vector<point>* data, vector<vertex>* tree, int step) {
	//make vertices
	(*tree)[step].data_start = start;
	(*tree)[step].data_end = end;
	for (int i = 0; i < (*data)[0].hidden.size(); ++i) {
		(*tree)[step].hidden_center.push_back(0); // new coordinate
		for (int j = start; j < end; ++j) {
			(*tree)[step].hidden_center[i] += (*data)[j].hidden[i]; //sum around start-end
		}
		(*tree)[step].hidden_center[i] /= end - start;//normalize
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

	// searching for the largest dimension-span {dim}
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

	for (int i = start; i < end; ++i) {
		(*data)[i].sorted_param = dim;
	}

	// medianth-stat partial sort and recursive call of the general fuction
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
// next functions are callable in dll
////////////////////////////////////////////////////////////

// Next fuction takes given signal and finds the nearest point in database saved by (build_tree) fuction.
// The found solution will be in array (double* nearest), also the fuction returns various info about considered clusters 
// (each is either a single point or a cluster that was discarded as a whole, i.e. they represent the highest level in the
// tree reached along each branch). Below we denote them "final clusters".
// The function works only with preliminary allocated arrays.
// number_points, hidden_size, measured_size - sizes of dataset, 
// number of object parameters (characteristics) and number of points in signal, respectively.
// measure - measured points in format {measured...}
// above array are preliminary constructed, following ones are resulting outputs
// centers - parameters, corresponding to centers of final clusters {hidden...}^K, K is discussed below 
// signal - signals, corresponding to centers of clusters probs {measured...}^K
// dist_s - distances to clusters (needed to calculate probs)
// weights - numbers of points in clusters, further used as weights of these clusters
// nearest - nearest found point
// ratio_perm - the maximum ratio of R and dist for rejecting cluster, otherwise we continue to split
// is_build_heap - whether we need the heap building (for optimal choice of the next cluster, normally used)
// max_clust_num - maximum number of returned clusters (for which memory is allocated). If larger number of clusters 
//   is actually considered during search, the error is returned.
// return - number of final clusters (K), satisfying K<=number_points and K<=max_clust_num
DllExport int c_probs(int number_points, int hidden_size, int measured_size,
	double* measure, double* centers, double* signal, double* dist_s, int* weights, double* nearest,
	double ratio_perm, int is_build_heap, bool load, int max_clust_num) {

	//load if needed (only from the default folder), normally not used
	if (load) {
		FILE* fout = fopen("built_tree_data.bin", "rb");
		FILE* fconf = fopen("common_info.conf", "rt");
		long number_points;
		fscanf(fconf, "%ld", &number_points);
		fscanf(fconf, "%d", &hidden_size);
		fscanf(fconf, "%d", &measured_size);

		m_size = measured_size;
		h_size = hidden_size;
		tr_size = 2 << (int(logb(number_points)) + 1);
		tree_starts_saved = new int[tr_size * 2];
		tree_saved = new double[(1 + hidden_size + measured_size) * tr_size];

		fread(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fread(tree_saved, sizeof(double), tr_size * (hidden_size + measured_size + 1), fout);
		fclose(fconf);
		fclose(fout);
		return 0;
	}

	// asserts
	if (hidden_size != h_size)
		return -1;
	if (measured_size != m_size)
		return -2;
	if (tree_saved == 0)
		return -3;

	// changing vector to standard array may potentially improve performance 
	vector<clust> cl_stack;// stack of clusters for testing
	clust first = { 0, count_dist(tree_saved + 1 + hidden_size, measure, measured_size) };
	cl_stack.push_back(first);
	if (is_build_heap)
		make_heap(cl_stack.begin(), cl_stack.end());

	int probs_counter = 0;
	double upper_bound_dist = DBL_MAX;
	double smallest_dist = DBL_MAX;
	int comps_num = 0;
	while (cl_stack.size() > 0) {
		comps_num++;
		int split = 1;
		clust cur;
		if (is_build_heap)
			cur = cl_stack[0];
		else
			cur = cl_stack[cl_stack.size() - 1];
		double dist = count_dist(tree_saved + (hidden_size + measured_size + 1) * cur.num + 1 + hidden_size,
			measure, measured_size);

		// Work with the chosen element
		//*****

		// do the following if cluster is to be discarded
		if ((upper_bound_dist < dist - (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0] &&
			ratio_perm*dist > (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0]) ||
			tree_starts_saved[2 * (cur.num) + 1] - tree_starts_saved[2 * (cur.num)] == 1) { // here test of splitting
			if (max_clust_num <= probs_counter)
				return -4;
			split = 0;
			dist_s[probs_counter] = dist*dist;
			weights[probs_counter] = tree_starts_saved[2 * cur.num + 1] - tree_starts_saved[2 * cur.num];
			copy(tree_saved + cur.num * (hidden_size + measured_size + 1) + 1,
				tree_saved + cur.num * (hidden_size + measured_size + 1) + 1 + hidden_size,
				centers + probs_counter * hidden_size);
			copy(tree_saved + cur.num * (hidden_size + measured_size + 1) + 1 + hidden_size,
				tree_saved + (cur.num + 1) * (hidden_size + measured_size + 1),
				signal + probs_counter * hidden_size);
			probs_counter++;
		}

		// new upper bound
		if (upper_bound_dist > dist + (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0])
			upper_bound_dist = dist + (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0];
		//*****
		if (smallest_dist > dist + *(tree_saved + cur.num * (hidden_size + measured_size + 1))) {
			smallest_dist = dist + *(tree_saved + cur.num * (hidden_size + measured_size + 1));
			copy(tree_saved + cur.num * (hidden_size + measured_size + 1) + 1,
				tree_saved + (cur.num + 1) * (hidden_size + measured_size + 1),
				nearest);
		}
		//delete it from stack
		if (is_build_heap)
			pop_heap(cl_stack.begin(), cl_stack.end());
		cl_stack.pop_back();
		//put new elements
		if (split) {
			clust sec = { 2 * cur.num + 1, count_dist(tree_saved + (2 * cur.num + 1) * (1 + hidden_size + measured_size)
				+ 1 + hidden_size, measure, measured_size) /*- tree[(2 * cur.num + 1)*(1 + hidden_size + measured_size)]*/ };
			clust third = { 2 * cur.num + 2, count_dist(tree_saved + (2 * cur.num + 2) * (1 + hidden_size + measured_size)
				+ 1 + hidden_size, measure, measured_size) /*- tree[(2 * cur.num + 2)*(1 + hidden_size + measured_size)]*/ };
			if (is_build_heap) {
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
}


//description is the same as for c_probs
//this version is just a shorter version of the above one, it doesn't return average clusters signals
DllExport int c_probs_fast(int number_points, int hidden_size, int measured_size,
	double* measure, double* centers, double* dist_s, int* weights, double* nearest,
	double ratio_perm, int is_build_heap, bool load, int max_clust_num) {

	//load if needed
	if (load) {
		FILE* fout = fopen("built_tree_data.bin", "rb");
		FILE* fconf = fopen("common_info.conf", "rt");
		long number_points;
		fscanf(fconf, "%ld", &number_points);
		fscanf(fconf, "%d", &hidden_size);
		fscanf(fconf, "%d", &measured_size);

		m_size = measured_size;
		h_size = hidden_size;
		tr_size = 2 << (int(logb(number_points)) + 1);
		tree_starts_saved = new int[tr_size * 2];
		tree_saved = new double[(1 + hidden_size + measured_size) * tr_size];

		fread(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fread(tree_saved, sizeof(double), tr_size * (hidden_size + measured_size + 1), fout);
		fclose(fconf);
		fclose(fout);
		return 0;
	}

	// asserts
	if (hidden_size != h_size)
		return -1;
	if (measured_size != m_size)
		return -2;
	if (tree_saved == 0)
		return -3;

	vector<clust> cl_stack;// stack of clusters for testing
	clust first = { 0, count_dist(tree_saved + 1 + hidden_size, measure, measured_size) };
	cl_stack.push_back(first);
	if (is_build_heap)
		make_heap(cl_stack.begin(), cl_stack.end());

	int probs_counter = 0;
	double upper_bound_dist = DBL_MAX;
	double smallest_dist = DBL_MAX;
	int comps_num = 0;
	while (cl_stack.size() > 0) {
		comps_num++;
		int split = 1;
		clust cur;
		if (is_build_heap)
			cur = cl_stack[0];
		else
			cur = cl_stack[cl_stack.size() - 1];
		double dist = count_dist(tree_saved + (hidden_size + measured_size + 1) * cur.num + 1 + hidden_size,
			measure, measured_size);

		// Work with chosen element
		//******

		// do the following if cluster is to be discarded
		if ((upper_bound_dist < dist - (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0] &&
			ratio_perm * dist >(tree_saved + (hidden_size + measured_size + 1) * cur.num)[0]) ||
			tree_starts_saved[2 * (cur.num) + 1] - tree_starts_saved[2 * (cur.num)] == 1) { // here test of splitting
			if (max_clust_num <= probs_counter)
				return -4;
			split = 0;
			dist_s[probs_counter] = dist*dist;
			weights[probs_counter] = tree_starts_saved[2 * cur.num + 1] - tree_starts_saved[2 * cur.num];
			//next also may be deleted for option wich return neither centers of thrown clusters
			copy(tree_saved + cur.num * (hidden_size + measured_size + 1) + 1,
				tree_saved + cur.num * (hidden_size + measured_size + 1) + 1 + hidden_size,
				centers + probs_counter * hidden_size);
			probs_counter++;
		}

		// new upper bound
		if (upper_bound_dist > dist + (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0])
			upper_bound_dist = dist + (tree_saved + (hidden_size + measured_size + 1) * cur.num)[0];
		//*****
		if (smallest_dist > dist + *(tree_saved + cur.num * (hidden_size + measured_size + 1))) {
			smallest_dist = dist + *(tree_saved + cur.num * (hidden_size + measured_size + 1));
			copy(tree_saved + cur.num * (hidden_size + measured_size + 1) + 1,
				tree_saved + (cur.num + 1) * (hidden_size + measured_size + 1),
				nearest);
		}
		//delete it from stack
		if (is_build_heap)
			pop_heap(cl_stack.begin(), cl_stack.end());
		cl_stack.pop_back();
		//put new elements
		if (split) {
			clust sec = { 2 * cur.num + 1, count_dist(tree_saved + (2 * cur.num + 1) * (1 + hidden_size + measured_size)
				+ 1 + hidden_size, measure, measured_size) /*- tree[(2 * cur.num + 1)*(1 + hidden_size + measured_size)]*/ };
			clust third = { 2 * cur.num + 2, count_dist(tree_saved + (2 * cur.num + 2) * (1 + hidden_size + measured_size)
				+ 1 + hidden_size, measure, measured_size) /*- tree[(2 * cur.num + 2)*(1 + hidden_size + measured_size)]*/ };
			if (is_build_heap) {
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
}

// Function below takes signals database and construct a pseudo-tree (of kd structure) 
// in the internal dll memory, it will be used by another fuction (c_probs or c_probs_fast).
// Also the fuction may be used to open the previouly constructed tree (set use_saved = 1) or save it (set save = 1).
// function works with allocated arrays
// data - input data in format {hidden..., measured...}^N
// number_points, hidden_size, measured_size - sizes of dataset, 
// number of object parameters (characteristics) and number of points in signal, respectively.
// (number_points can be smaller than total number of rows in data, only the former will be used)
// save - bool, whether you wish to save your tree as file
// use_saved - bool, whether you want to load a tree from folder instead of building it
// folder_name - folder to use for saving or loading the tree
DllExport int32_t build_tree(double* data, int32_t number_points,
	int32_t hidden_size, int32_t measured_size, int32_t save, 
	int32_t use_saved, char* folder_name) {

	//construct saving(reading) tree filename
	std::string f_n1(folder_name);
	f_n1 += "\\built_tree_data.bin";
	std::string f_n2(folder_name);
	f_n2 += "\\common_info.conf";

	if (use_saved) {
		FILE* fout = fopen(f_n1.c_str(), "rb");
		FILE* fconf = fopen(f_n2.c_str(), "rt");
		if ((fout == 0) || (fconf == 0))
			return -5;
		fscanf(fconf, "%ld", &number_points);
		fscanf(fconf, "%d", &hidden_size);
		fscanf(fconf, "%d", &measured_size);

		m_size = measured_size;
		h_size = hidden_size;
		tr_size = 2 << (int(logb(number_points)) + 1);
		tree_starts_saved = new int[tr_size * 2];
		tree_saved = new double[(1 + hidden_size + measured_size) * tr_size];

		fread(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fread(tree_saved, sizeof(double), tr_size * (hidden_size + measured_size + 1), fout);
		fclose(fconf);
		fclose(fout);
		return 1;
	}

	//memory needed to allocate (with some additional capacity)
	tr_size = 2 << (int(logb(number_points)) + 1);

	vector<point>* data_vect = transform_measured_data(data, number_points, hidden_size, measured_size);
	vector<vertex>* tree_vect = new vector<vertex>(tr_size);

	// only this function call build tree itself
	// all other in this fuction is just saving or preparation
	build(0, number_points, data_vect, tree_vect, 0);

	tree_starts_saved = new int[tr_size * 2];
	tree_saved = new double[(1 + hidden_size + measured_size) * tr_size];

	m_size = measured_size;
	h_size = hidden_size;

	for (int i = 0; i < tree_vect->size(); ++i) {
		tree_starts_saved[2 * i] = (*tree_vect)[i].data_start;
		tree_starts_saved[2 * i + 1] = (*tree_vect)[i].data_end;

		tree_saved[i * (1 + hidden_size + measured_size)] = (*tree_vect)[i].R;
		for (int j = 0; j < (*tree_vect)[i].hidden_center.size(); ++j) {
			tree_saved[i * (1 + hidden_size + measured_size) + 1 + j] = (*tree_vect)[i].hidden_center[j];
		}
		for (int j = 0; j < (*tree_vect)[i].measured_center.size(); ++j) {
			tree_saved[i * (1 + hidden_size + measured_size) + 1 + hidden_size + j] = (*tree_vect)[i].measured_center[j];
		}
	}

	if (save) {
		FILE* fout = fopen(f_n1.c_str(), "wb");
		FILE* fconf = fopen(f_n2.c_str(), "wt");
		fprintf(fconf, "%ld ", number_points);
		fprintf(fconf, "%d ", hidden_size);
		fprintf(fconf, "%d ", measured_size);
		fwrite(tree_starts_saved, sizeof(int), tr_size * 2, fout);
		fwrite(tree_saved, sizeof(double), tr_size * (hidden_size + measured_size + 1), fout);
		fclose(fconf);
		fclose(fout);	
	}

	return 1;
}


//////////////////////////////////////////////////
//in the following there are only tests, which are not used in the dll aplication (only for debug)
////////////////////////////////////////////////////

// Fuction for test of saving tree data
DllExport int get_msize(int var) {
	return m_size;
}

double* load_binary_data(const char* name, int* num, int* hidden, int* measured) {
	FILE* f = fopen(name, "rb");
	fread(num, sizeof(int), 1, f);
	fread(hidden, sizeof(int), 1, f);
	fread(measured, sizeof(int), 1, f);
	int double_nums = *num * (*hidden + *measured);
	double* data = new double[double_nums];
	auto r = fread(data, sizeof(double), double_nums, f);
	return data;
}

void build_and_save_tree() {
	int hid_size = 4;
	int mes_size = 61;
	int database_size = 300000;

	// read full database, normed not in another program
	double* data = load_binary_data("databasebytes_new", &database_size, &hid_size, &mes_size);

	int upperbound_tree = 2 << (int(logb(database_size)) + 1);

	//build_tree(data, database_size, hid_size, mes_size, true, 0);

	return;
}
