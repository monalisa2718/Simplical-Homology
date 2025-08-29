#include<iostream>
#include<string>
#include<vector>
#include<cassert>
#include<algorithm>
#include<fstream>
#include <map>
#include <set>
#include <queue>

using namespace std;

vector<vector<int>> mat_inv;
class generator_of_module {
public:
	vector<int> name;
	int coefficient = 0;
};
class simplicial_chain {
public:
	vector<vector<int>> simplical_complex;
	vector<vector<int>> simplical_precomplex;
	vector<int> rank;
	vector<vector<vector<int>>> base;
	vector<vector<int>> boundary_matrix;
	vector<vector<int>> coboundary_matrix;
	void set_rank_base() {
		int max_degree = 0;
		rank.resize(0);
		for (int i = 0; i < simplical_complex.size(); i++) {
			if (max_degree < simplical_complex[i].size()) {
				max_degree = simplical_complex[i].size();
			}
		}
		for (int i = 1; i <= max_degree; i++) {
			vector<vector<int>> vvs;
			for (int j = 0; j < simplical_complex.size(); j++) {
				if (simplical_complex[j].size() == i) {
					vvs.push_back(simplical_complex[j]);
				}
			}
			base.push_back(vvs);
			rank.push_back(vvs.size());
		}
	};
	void show_base() {
		for (int i = 0; i < base.size(); i++) {
			for (int j = 0; j < base[i].size(); j++) {
				cout << "{";
				for (int k = 0; k < base[i][j].size(); k++) {
					cout << base[i][j][k];
					if (k != base[i][j].size() - 1) cout << ",";
				}
				cout << "}" << j << ", ";
				if ((j + 1) % 10 == 0) cout << endl;
			}
			cout << endl << endl;
		}
	}
	void show_base_text(ofstream* pttext) {
		for (int i = 0; i < base.size(); i++) {
			for (int j = 0; j < base[i].size(); j++) {
				*pttext << "{";
				for (int k = 0; k < base[i][j].size(); k++) {
					*pttext << base[i][j][k];
					if (k != base[i][j].size() - 1) *pttext << ",";
				}
				*pttext << "}" << j << ", ";
				if ((j + 1) % 10 == 0) *pttext << endl;
			}
			*pttext << endl << endl;
		}
	}
	void show_rank() {
		for (int i = 0; i < rank.size(); i++) {
			cout << "Chain's rank of degree " << i << " is " << rank[i] << "." << endl;
		}
	}
	void show_coboundary_matrix() {
		for (vector<int> vi : coboundary_matrix) {
			for (int i : vi) {
				cout << i << " ";
			}
			cout << endl;
		}
	}
	void show_boundary_matrix() {
		for (vector<int> vi : boundary_matrix) {
			for (int i : vi) {
				cout << i << " ";
			}
			cout << endl;
		}
	}
	void set_coboundary_matrix(int n) {
		//この関数でn-1次のコチェインからn次のコチェインへのコバウンダリー写像を作ります。
		set_rank_base();
		coboundary_matrix.resize(0);
		coboundary_matrix.resize(rank[n]);
		for (int i = 0; i < coboundary_matrix.size(); i++) {
			coboundary_matrix[i].resize(rank[n - 1], 0);
		}
		for (int j = 0; j < rank[n]; j++) {
			for (generator_of_module gm : boundary(base[n][j])) {
				for (int i = 0; i < rank[n - 1]; i++) {
					if (base[n - 1][i] == gm.name) {
						coboundary_matrix[j][i] = gm.coefficient;
						break;
					}
				}
			}
		}
	}
	void set_boundary_matrix(int n) {
		set_coboundary_matrix(n);
		vector<int> vi;
		vi.resize(coboundary_matrix.size());
		boundary_matrix.resize(0);
		boundary_matrix.resize(coboundary_matrix[0].size(), vi);
		for (int i = 0; i < coboundary_matrix.size(); i++) {
			for (int j = 0; j < coboundary_matrix[0].size(); j++) {
				boundary_matrix[j][i] = coboundary_matrix[i][j];
			}
		}
	}
	void complexed() {
		simplical_complex.resize(0);
		for (vector<int> vs : simplical_precomplex) {
			for (int bit = 0; bit < (1 << vs.size()); bit++) {
				vector<int> vs_a;
				for (int i = 0; i < vs.size(); i++) {
					if (bit & (1 << i)) {
						vs_a.push_back(vs[i]);
					}
				}
				if (find(simplical_complex.begin(), simplical_complex.end(), vs_a) == simplical_complex.end()) {
					simplical_complex.push_back(vs_a);
				}
			}
		}
	}//not complited
private:
	vector<generator_of_module> boundary(vector<int> element_of_chain) {
		vector<generator_of_module> vg;
		element_of_chain.size();
		for (int i = 0; i < element_of_chain.size(); i++) {
			vg.resize(vg.size() + 1);
			for (int j = 0; j < element_of_chain.size(); j++) {
				if (i != j) {
					vg[i].name.push_back(element_of_chain[j]);
				}
			}
			vg[i].coefficient = pow(-1, i);
		}
		return vg;
	}
};

void show_matrix(vector<vector<int>> matrix) {
	cout << "{";
	for (int i = 0; i < matrix.size(); i++) {
		cout << "{";
		for (int j = 0; j < matrix[0].size(); j++) {
			cout << matrix[i][j];
			if (j != (matrix[0].size() - 1)) {
				cout << ",";
			}
		}
		cout << "}";
		if (i != matrix.size() - 1) {
			cout << "," << endl;
		}
	}
	cout << "}" << endl << endl;
}
void show_matrix_text(vector<vector<int>> matrix, ofstream* pttext) {
	cout << "{";
	for (int i = 0; i < matrix.size(); i++) {
		*pttext << "{";
		for (int j = 0; j < matrix[0].size(); j++) {
			*pttext << matrix[i][j];
			if (j != (matrix[0].size() - 1)) {
				*pttext << ",";
			}
		}
		*pttext << "}";
		if (i != matrix.size() - 1) {
			*pttext << "," << endl;
		}
	}
	*pttext << "}" << endl << endl;
}
void show_matrix_sage(vector<vector<int>> matrix) {
	cout << "[";
	for (int i = 0; i < matrix.size(); i++) {
		cout << "[";
		for (int j = 0; j < matrix[0].size(); j++) {
			cout << matrix[i][j];
			if (j != (matrix[0].size() - 1)) {
				cout << ",";
			}
		}
		cout << "]";
		if (i != matrix.size() - 1) {
			cout << ",";
		}
	}
	cout << "]" << endl;
}
void slide_push_back(vector<vector<int>>* ptvvi, int num) {
	const int n = (*ptvvi).size();
	const int m = (*ptvvi)[n - 1].size();
	(*ptvvi).resize(n + 1);
	for (int i = 1; i < m; i++) {
		(*ptvvi)[n].push_back((*ptvvi)[n - 1][i]);
	}
	(*ptvvi)[n].push_back((*ptvvi)[n][0] + num);
}
void show_vvi(vector<vector<int>> vvi) {
	for (int i = 0; i < vvi.size(); i++) {
		for (int j = 0; j < vvi[i].size(); j++) {
			cout << vvi[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
void sort_vvi(vector<vector<int>>* ptvvi) {
	for (vector<int>& vi : *ptvvi) {
		sort(vi.begin(), vi.end());
	}
}
void show_vi(vector<int> vi) {
	const int n = vi.size();
	cout << "{";
	for (int i = 0; i < n; i++) {
		cout << vi[i] ;
		if (i != n - 1) cout << ",";
	}
	cout << "}" << endl;
}
void show_orient(vector<generator_of_module> vg) {
	for (generator_of_module gm : vg) {
		cout << gm.coefficient << " : {";
		for (int i = 0; i < gm.name.size(); i++) {
			cout << gm.name[i] ;
			if (i != gm.name.size() - 1) {
				cout << ",";
			}
		}
		cout << "}" << endl;
	}
}
vector<vector<int>> times_circle_simplex(vector<int> vi, int num) {
	vector<vector<int>> vvi;
	const int n = vi.size();
	vvi.push_back(vi);
	vvi[0].push_back(vvi[0][0] + num);
	for (int i = 1; i < n; i++) {
		slide_push_back(&vvi, num);
	}
	for (int i = 0; i < n; i++) {
		slide_push_back(&vvi, -num);
	}
	return vvi;
}
vector<vector<int>> times_circle_three_simplex(vector<int> vi, int num, int grade) {
	vector<vector<int>> vvi;
	const int n = vi.size();
	vvi.push_back(vi);
	vvi[0].push_back(vvi[0][0] + num);
	for (int i = 1; i < grade * n; i++) {
		slide_push_back(&vvi, num);
	}
	for (int i = 0; i < n; i++) {
		slide_push_back(&vvi, - grade * num);
	}
	return vvi;
}
vector<vector<int>> vec_add(vector<vector<int>> vvi_a, vector<vector<int>>vvi_b) {
	for (auto vi : vvi_b) {
		vvi_a.push_back(vi);
	}
	return vvi_a;
}
vector<vector<int>> times_circle(vector<vector<int>> vvi) {
	vector<vector<int>> vvi_a;
	int m = (vvi)[0][0];
	int l = (vvi)[0][0];
	const int n = vvi.size();
	for (vector<int> vi : vvi) {
		for (int i : vi) {
			m = max(i, m);
			l = min(i, l);
		}
	}
	for (vector<int> vi : vvi) {
		vvi_a = vec_add(vvi_a, times_circle_simplex(vi, m - l + 1));
	}
	return vvi_a;
}
void times_circle_v(vector<vector<int>>* ptvvi) {
	vector<vector<int>> vvi_a;
	int m = (*ptvvi)[0][0];
	int l = (*ptvvi)[0][0];
	const int n = ptvvi-> size();
	for (vector<int> vi : *ptvvi) {
		for (int i : vi) {
			m = max(i, m);
			l = min(i, l);
		}
	}
	for (vector<int> vi : *ptvvi) {
		vvi_a = vec_add(vvi_a, times_circle_simplex(vi, m - l +1));
	}
	*ptvvi = vvi_a;
}
void times_circle_three_v(vector<vector<int>>* ptvvi, int grade) {
	vector<vector<int>> vvi_a;
	int m = (*ptvvi)[0][0];
	int l = (*ptvvi)[0][0];
	const int n = ptvvi->size();
	for (vector<int> vi : *ptvvi) {
		for (int i : vi) {
			m = max(i, m);
			l = min(i, l);
		}
	}
	for (vector<int> vi : *ptvvi) {
		vvi_a = vec_add(vvi_a, times_circle_three_simplex(vi, m - l + 1, grade));
	}
	*ptvvi = vvi_a;
}
vector<vector<int>> id(int n) {
	vector<vector<int>> vvi;
	vector<int> vi;
	vi.resize(n, 0);
	vvi.resize(n, vi);
	for (int i = 0; i < n; i++) {
		vvi[i][i] = 1;
	}
	return vvi;
}
vector<vector<int>> matrix_product(vector<vector<int>> matrix_a, vector<vector<int>> matrix_b) {
	if (matrix_a[0].size() != matrix_b.size()) {
		cout << "marix_product makes a mistake." << endl;
	}
	vector<vector<int>> vvi;
	vector<int> vi;
	vi.resize(matrix_b[0].size(), 0);
	vvi.resize(matrix_a.size(), vi);
	for (int i = 0; i < vvi.size(); i++) {
		for (int j = 0; j < vvi[0].size(); j++) {
			for (int k = 0; k < matrix_a[0].size(); k++) {
				vvi[i][j] += matrix_a[i][k] * matrix_b[k][j];
			}
		}
	}
	return vvi;
}
vector<vector<int>> matrix_product(vector<vector<int>> mat_a, vector<vector<int>> mat_b, vector<vector<int>> mat_c) {
	return matrix_product(mat_a, matrix_product(mat_b, mat_c));
}
vector<vector<int>> mat_scalar_prod(vector<vector<int>> vvi, int n) {
	vector<vector<int>> vvi_a(vvi.size());
	for (int i = 0; i < vvi_a.size(); i++) {
		for (int j : vvi[i]) {
			vvi_a[i].push_back(n * j);
		}
	}
	return vvi_a;
}
vector<vector<int>> mat_add(vector<vector<int>> vvi_a, vector<vector<int>> vvi_b) {
	vector<vector<int>> vvi_c;
	if (vvi_a.empty() || vvi_b.empty()) {
		cout << "mat_add error. empty matrix." << endl;
		return {};
	}
	if (vvi_a.size() != vvi_b.size() || vvi_a[0].size() != vvi_b[0].size()) {
		cout << "mat_add error. size don't match." << endl;
		return {};
	}
	assert(vvi_a.size() == vvi_b.size() && vvi_a[0].size() == vvi_b[0].size());
	vvi_c.resize(vvi_a.size());
	for (int i = 0; i < vvi_a.size(); i++) {
		for (int j = 0; j < vvi_a[0].size(); j++) {
			vvi_c[i].push_back(vvi_a[i][j] + vvi_b[i][j]);
		}
	}
	return vvi_c;
}
vector<vector<int>> add_column(vector<vector<int>>* matrix, int n, int m, int l) {//n列目をm列目にl倍して加える関数。
	vector<vector<int>> vvi = id((*matrix)[0].size());
	vvi[n][m] += l;
	for (int i = 0; i < (*matrix).size(); i++) {
		(*matrix)[i][m] += ((*matrix)[i][n] * l);
	}
	return vvi;
}
void add_column_v(vector<vector<int>>* matrix, int n, int m, int l) {//n列目をm列目にl倍して加える関数。
	vector<vector<int>> vvi = id((*matrix)[0].size());
	vvi[n][m] += l;
	for (int i = 0; i < (*matrix).size(); i++) {
		(*matrix)[i][m] += ((*matrix)[i][n] * l);
	}
}
vector<vector<int>> add_row(vector<vector<int>>* matrix, int n, int m, int l) {//n行目をm行目にl倍して加える関数。
	vector<vector<int>> vvi = id((*matrix).size());
	vvi[m][n] += l;
	for (int i = 0; i < (*matrix)[0].size(); i++) {
		(*matrix)[m][i] = (*matrix)[n][i] * l + (*matrix)[m][i];
	}
	return vvi;
}
void add_row_v(vector<vector<int>>* matrix, int n, int m, int l) {//n行目をm行目にl倍して加える関数。
	vector<vector<int>> vvi = id((*matrix).size());
	vvi[m][n] += l;
	for (int i = 0; i < (*matrix)[0].size(); i++) {
		(*matrix)[m][i] = (*matrix)[n][i] * l + (*matrix)[m][i];
	}
}
vector<vector<int>> exchange_column(vector<vector<int>>* ptmatrix, int n, int m) {
	vector<vector<int>> vvi = id((*ptmatrix)[0].size());
	vvi[n][n] = 0;
	vvi[m][m] = 0;
	vvi[m][n] = 1;
	vvi[n][m] = 1;
	int j;
	for (int i = 0; i < (*ptmatrix).size(); i++) {
		j = (*ptmatrix)[i][n];
		(*ptmatrix)[i][n] = (*ptmatrix)[i][m];
		(*ptmatrix)[i][m] = j;
	}
	return vvi;
}
void exchange_column_v(vector<vector<int>>* ptmatrix, int n, int m) {
	vector<vector<int>> vvi = id((*ptmatrix)[0].size());
	vvi[n][n] = 0;
	vvi[m][m] = 0;
	vvi[m][n] = 1;
	vvi[n][m] = 1;
	int j;
	for (int i = 0; i < (*ptmatrix).size(); i++) {
		j = (*ptmatrix)[i][n];
		(*ptmatrix)[i][n] = (*ptmatrix)[i][m];
		(*ptmatrix)[i][m] = j;
	}
}
vector<vector<int>> exchange_row(vector<vector<int>>* ptmatrix, int n, int m) {
	vector<vector<int>> vvi = id((*ptmatrix).size());
	vvi[n][n] = 0;
	vvi[m][m] = 0;
	vvi[m][n] = 1;
	vvi[n][m] = 1;
	int j;
	for (int i = 0; i < (*ptmatrix)[0].size(); i++) {
		j = (*ptmatrix)[n][i];
		(*ptmatrix)[n][i] = (*ptmatrix)[m][i];
		(*ptmatrix)[m][i] = j;
	}
	return vvi;
}
void exchange_row_v(vector<vector<int>>* ptmatrix, int n, int m) {
	vector<vector<int>> vvi = id((*ptmatrix).size());
	vvi[n][n] = 0;
	vvi[m][m] = 0;
	vvi[m][n] = 1;
	vvi[n][m] = 1;
	int j;
	for (int i = 0; i < (*ptmatrix)[0].size(); i++) {
		j = (*ptmatrix)[n][i];
		(*ptmatrix)[n][i] = (*ptmatrix)[m][i];
		(*ptmatrix)[m][i] = j;
	}
}
vector<vector<int>> transpose(vector<vector<int>> mat) {
	vector<vector<int>> vvi;
	vvi.resize(mat[0].size());
	for (int i = 0; i < vvi.size(); i++) {
		vvi[i].resize(mat.size(), 0);
		for (int j = 0; j < mat.size(); j++) {
			vvi[i][j] = mat[j][i];
		}
	}
	return vvi;
}
vector<vector<vector<int>>> Smith_normalize(vector<vector<int>> matrix) {
	//この関数はバウンダリーやコバウンダリー等の特殊な行列でしか動きません。
	vector<vector<vector<int>>> vvvi;
	vvvi.resize(3);
	vvvi[0] = matrix;
	vvvi[1] = id(matrix.size());
	vvvi[2] = id(matrix[0].size());

	for (int search_number = 0; (search_number < vvvi[0].size()) && (search_number < vvvi[0][0].size()); search_number++) {
		int find_one[3] = { 0,0,0 };
		for (int i = search_number; i < matrix.size(); i++) {
			find_one[2] = 0;
			for (int j = search_number; j < matrix[0].size(); j++) {
				if (abs(vvvi[0][i][j]) == 1) {
					find_one[0] = i;
					find_one[1] = j;
					find_one[2]++;
					if (vvvi[0][i][j] == -1) {
						for (int k = 0; k < vvvi[0][0].size(); k++) {
							vvvi[0][i][k] = -1 * vvvi[0][i][k];
						}
						for (int k = 0; k < vvvi[1][0].size(); k++) {
							vvvi[1][i][k] = -1 * vvvi[1][i][k];
						}
					}
					break;
				}
			}
			if (find_one[2] != 0)break;
		}
		if (find_one[2] != 0) {
			//vvvi[1] = matrix_product(exchange_row(&(vvvi[0]), search_number, find_one[0]), vvvi[1]);
			exchange_row_v(&(vvvi[0]), search_number, find_one[0]);
			exchange_row_v(&(vvvi[1]), search_number, find_one[0]);
			//vvvi[2] = matrix_product(vvvi[2], exchange_column(&(vvvi[0]), search_number, find_one[1]));
			exchange_column_v(&(vvvi[0]), search_number, find_one[1]);
			exchange_column_v(&(vvvi[2]), search_number, find_one[1]);
			for (int i = search_number + 1; i < vvvi[0].size(); i++) {
				//vvvi[1] = matrix_product(add_row(&(vvvi[0]), search_number, i, -1 * vvvi[0][i][search_number]), vvvi[1]);
				add_row_v(&(vvvi[1]), search_number, i, -1 * vvvi[0][i][search_number]);
				add_row_v(&(vvvi[0]), search_number, i, -1 * vvvi[0][i][search_number]);
			}
			for (int i = search_number + 1; i < vvvi[0][0].size(); i++) {
				//vvvi[2] = matrix_product(vvvi[2], add_column(&(vvvi[0]), search_number, i, -1 * vvvi[0][search_number][i]));
				add_column_v(&(vvvi[2]), search_number, i, -1 * vvvi[0][search_number][i]);
				add_column_v(&(vvvi[0]), search_number, i, -1 * vvvi[0][search_number][i]);
			}
		}
	}
	return vvvi;
}
vector<vector<vector<int>>> Smith_normalize_L(vector<vector<int>> mat) {
	vector<vector<vector<int>>> vvvi;
	vvvi.push_back(mat);
	vvvi.push_back(id(mat.size()));
	show_matrix(vvvi[0]);
	show_matrix(vvvi[1]);
	int find_one[3] = { 0,0,2 };
	for (int i = 0; i < mat[0].size() && i < mat.size(); i++) {
		for (int j = i; j < mat.size(); j++) {
			find_one[0] = 0;
			find_one[1] = 0;
			find_one[2] = 2;
			if (vvvi[0][j][i] == 1 || vvvi[0][j][i] == -1) {
				find_one[0] = j;
				find_one[1] = i;
				find_one[2] = vvvi[0][j][i];
				if (vvvi[0][j][i] == -1) {
					for (int k = 0; k < mat[0].size(); k++) {
						vvvi[0][j][k] = -1 * vvvi[0][j][k];
						vvvi[1][j][k] = -1 * vvvi[1][j][k];
					}
				}
				exchange_row_v(&(vvvi[1]), i, find_one[0]);
				exchange_row_v(&(vvvi[0]), i, find_one[0]);
				break;
			}
			if (find_one[2] == 2 && j == mat.size() - 1) {
				cout << "no one" << endl;
			}
		}
		for (int j = 0; j < mat.size() && find_one[2] != 2; j++) {
			if (j != i) {
				add_row_v(&(vvvi[1]), i, j, -1 * vvvi[0][j][i]);
				add_row_v(&(vvvi[0]), i, j, -1 * vvvi[0][j][i]);
			}
		}
		cout << "check" << endl;
		show_matrix(vvvi[0]);
		show_matrix(vvvi[1]);
		show_matrix(matrix_product(vvvi[1], mat));
	}
	return vvvi;
}
vector<vector<vector<int>>> inverse(vector<vector<int>> matrix) {
	//この関数はバウンダリーやコバウンダリー等の特殊な行列でしか動きません。
	if (true)
	{//implement by Ryo Takahashi
		vector<vector<vector<int>>>vvvi(2);
		if (matrix.empty())
		{
			cout << "Empty matrix" << endl;
			return {};//return empty
		}
		if (matrix.size() != matrix[0].size())
		{
			cout << "Not a square matrix" << endl;
			return {};
		}
		assert(matrix.size() == matrix[0].size());
		const int n = matrix.size();//size of the square matrix
		vvvi[0] = matrix;
		vvvi[1] = id(n);
		for (int s = 0; s < n; s++)
		{
			int i = s;
			while (i < n && abs(vvvi[0][i][s]) != 1)i++;
			if (i == n)
			{
				cout << "no one at s =" << s << endl;
				for (int x = 0; x < n; x++)
				{
					for (int y = 0; y < n; y++)cout << vvvi[0][x][y] << (y + 1 == n ? "\n" : " ");
				}
			}
			if (i > s) {
				swap(vvvi[0][i], vvvi[0][s]);
				swap(vvvi[1][i], vvvi[1][s]);
			}
			if (vvvi[0][s][s] == -1)
			{
				for (int j = 0; j < n; j++)
				{
					vvvi[0][s][j] = -vvvi[0][s][j];
					vvvi[1][s][j] = -vvvi[1][s][j];
				}
			}
			assert(vvvi[0][s][s] == 1);
			for (int i = 0; i < n; i++)if (i != s && vvvi[0][i][s] != 0)
			{
				const int coef = -vvvi[0][i][s];
				for (int j = 0; j < n; j++)
				{
					vvvi[0][i][j] += vvvi[0][s][j] * coef;
					vvvi[1][i][j] += vvvi[1][s][j] * coef;
				}
			}
		}
		return vvvi;
	}
	if (false) {
		vector<vector<vector<int>>> vvvi;
		vvvi.resize(2);
		vvvi[0] = matrix;
		vvvi[1] = id(matrix.size());
		for (int search_number = 0; search_number < vvvi[0].size(); search_number++) {
			if (matrix.size() != matrix[0].size()) {
				cout << "inverse can't be used for not square matricis." << endl;
				break;
			}
			cout << search_number << endl;
			int find_one[3] = { 0,0,0 };
			for (int i = search_number; i < matrix.size(); i++) {
				find_one[2] = 0;
				if (vvvi[0][i][search_number] == 1 || vvvi[0][i][search_number] == -1) {
					find_one[0] = i;
					find_one[1] = search_number;
					find_one[2]++;
					break;
				}//-1倍していない
				if (find_one[2] == 0 && i == matrix.size() - 1) {
					cout << "one can't be finded." << endl;
				}
			}
			if (find_one[2] != 0) {
				vvvi[1] = matrix_product(exchange_row(&(vvvi[0]), search_number, find_one[0]), vvvi[1]);
				for (int i = 0; i < vvvi[0].size(); i++) {
					if (i != search_number) {
						//vvvi[1] = matrix_product(add_row(&(vvvi[0]), search_number, i, -1 * vvvi[0][find_one[0]][search_number] * vvvi[0][i][search_number]), vvvi[1]);
						vvvi[1] = matrix_product(add_row(&(vvvi[0]), search_number, i, -1 * vvvi[0][i][search_number]), vvvi[1]);
					}
				}
			}
		}
		return vvvi;
	}
}
vector<vector<int>> same_column(vector<vector<int>> mat_A, vector<vector<int>> mat_B) {
	vector<vector<int>> vvi;
	vector<vector<int>> vvi_a;
	vector<vector<int>> vvi_b;
	if (mat_A.size() != mat_B.size()) {
		cout << "same_column can't work correctly." << endl;
		return vvi;
	}
	vvi_a = transpose(mat_A);
	vvi_b = transpose(mat_B);
	for (vector<int> vi_A : vvi_a) {
		for (vector<int> vi_B : vvi_b) {
			if (vi_A == vi_B) {
				vvi.push_back(vi_A);
			}
		}
	}
	return vvi;
}
vector<generator_of_module> orient(vector<vector<int>> vvi) {
	if (vvi.empty()) {
		cout << "orient is givin empty vector." << endl;
		return {};
	}
	else {
		const int n = vvi[0].size();
		for (vector<int> vi : vvi) {
			if (vi.size() != n) {
				cout << "orient error. size is not appropriate." << endl;
				return {};
			}
		}
	}
	vector<generator_of_module> vg;
	vg.resize(vvi.size());
	vector<vector<int>> vvi_a;
	vector<vector<int>> vvi_b;
	const int n = vvi.size();
	const int l = vvi[0].size();
	int s;
	for (int i = 0; i < n; i++) {
		vg[i].name = vvi[i];
		vg[i].coefficient = 0;
	}
	vg[0].coefficient = 1;
	cout << "The orient of {";
	show_vi(vvi[0]);
	cout << "} is 1." << endl;
	assert(vvi.size() > 0);
	while (true) {
		s = 1;
		for (generator_of_module gm : vg) {
			s = s * gm.coefficient;
		}
		if (s != 0) break;
		for (int orient_ref = 0; orient_ref < n; orient_ref++) {
			if (vg[orient_ref].coefficient == 0) continue;
			vvi_a.resize(0);
			vvi_a.resize(l);
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < l; j++) {
					if (j != i) {
						vvi_a[i].push_back(vvi[orient_ref][j]);
					}
				}
			}
			for (int i = 0; i < n; i++) {
				if (i == orient_ref) continue;
				if (vg[i].coefficient != 0) continue;
				vvi_b.resize(0);
				vvi_b.resize(l);
				for (int j = 0; j < l; j++) {
					for (int k = 0; k < l; k++) {
						if (j != k) {
							vvi_b[j].push_back(vvi[i][k]);
						}
					}
				}
				s = 0;
				for (int j = 0; j < l; j++) {
					for (int k = 0; k < l; k++) {
						if (vvi_a[j] == vvi_b[k]) {
							vg[i].coefficient = vg[orient_ref].coefficient * (-1 + 2 * abs((j - k) % 2));
							s++;
							break;
						}
					}
					if (s == 1) break;
				}
			}
		}
	}
	return vg;
}

int main() {
	ofstream text("Homology_text");
	simplicial_chain ch;
	//ch.simplical_precomplex = { {0,1,2,3},{0,1,4,5},{0,2,4,6},{0,3,5,6},{1,2,5,7},{2,3,6,7},{1,4,7,8},{3,4,5,8},{2,4,6,9},{3,5,6,9},{1,5,7,9},{2,6,7,9},{4,5,8,9},{6,7,8,9},{1,3,7,8},{0,7,8,9} };
	//ch.simplical_precomplex = { {0,1,2,3,4},{1,2,3,4,5},{2,3,4,5,6},{3,4,5,6,7},{4,5,6,7,8},{5,6,7,8,9},{6,7,8,9,10},{7,8,9,10,11},{0,1,9,10,11},{0,1,2,10,11},{0,1,2,3,11} };
	ch.simplical_precomplex = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
	//ch.simplical_precomplex = { {0,1},{1,2},{0,2} };
	//times_circle_three_v(&(ch.simplical_precomplex),2);
	times_circle_three_v(&(ch.simplical_precomplex),2);
	times_circle_three_v(&(ch.simplical_precomplex),2);
	cout << endl;
	sort_vvi(&ch.simplical_precomplex);
	show_matrix(ch.simplical_precomplex);
	show_matrix_text(ch.simplical_precomplex, &text);
	ch.complexed();
	ch.set_rank_base();
	//ch.show_base();
	//ch.show_base_text(&text);
	vector<vector<vector<int>>> coboundary_one;
	vector<vector<vector<int>>> coboundary_two;
	vector<vector<vector<int>>> vvvi;
	vector<vector<int>> vvi_im;
	vector<vector<int>> vvi_ker;
	vector<vector<int>> A;
	vector<vector<int>> representative;
	vector<vector<int>> all_representative;
	vector<vector<int>> vvi_a;
	vector<vector<int>> vvi_b;
	vector<vector<int>> vvi_c;
	vector<int> vi;
	vector<generator_of_module> vg_orinent;
	int n;
	int m;
	int rank_ker;
	int rank_im = 0;
	int sum;
	ch.show_rank();
	vg_orinent = orient(ch.base[4]);
	show_orient(vg_orinent);
	if (true) {
		cout << "delta 1" << endl;
		text << "delta 1" << endl;
		ch.set_coboundary_matrix(2);
		coboundary_one = Smith_normalize(ch.coboundary_matrix);
		n = coboundary_one[1].size();
		vvi_a = inverse(coboundary_one[1])[1];
		vvi_im.resize(n);
		for (int i = 0; coboundary_one[0][i][i] == 1 ; i++) {
			for (int j = 0; j < n; j++) {
				vvi_im[j].push_back(vvi_a[j][i]);
			}
			rank_im++;
			if (i == min(coboundary_one[0].size(), coboundary_one[0][0].size()) - 1) {
				break;
			};
		}
		//show_matrix(ch.coboundary_matrix);
	}
	if (true) {
		cout << endl << "delta 2" << endl;
		text << endl << "delta 2" << endl;
		ch.set_coboundary_matrix(3);
		coboundary_two = Smith_normalize(ch.coboundary_matrix);
		vvi_ker.resize(coboundary_two[2].size());
		rank_ker = coboundary_two[2].size();
		cout << "rank_ker == " << rank_ker << endl;
		text << "rank_ker == " << rank_ker << endl;
		//show_matrix(coboundary_two[0]);
		cout << "min == " << min(coboundary_two[0].size(), coboundary_two[0][0].size()) << endl;
		for (int i = 0; coboundary_two[0][i][i] == 1; i++) {
			rank_ker--;
			if (i == min(coboundary_two[0].size(), coboundary_two[0][0].size()) - 1) {
				break;
			};
		}
		n = coboundary_two[2].size();
		for (int i = 0; i < rank_ker; i++) {
			for (int j = 0; j < n; j++) {
				vvi_ker[j].push_back(coboundary_two[2][j][n - i - 1]);
			}
		}
		//show_matrix(ch.coboundary_matrix);
	}
	if (true) {
		vvvi = Smith_normalize(vvi_ker);
		vvi_ker = matrix_product(vvi_ker, vvvi[2]);
		//vvi_im = matrix_product(vvi_im, vvvi[2]);
		//cout << "kar mat" << endl;
		//show_matrix(vvi_ker);
		vvi_a = matrix_product(vvvi[1], vvi_ker);
		//cout << "normalized kar mat" << endl;
		//show_matrix(vvi_a);
		vvi_b = matrix_product(vvvi[1], vvi_im);
		//cout << "normalized im mat" << endl;
		//show_matrix(vvi_b);
	}
	if (true) {
		A = vvi_b;
		A.resize(rank_ker);
		//cout << "A == " << endl;
		vvvi = Smith_normalize(A);
		vvi_im = matrix_product(vvi_im, vvvi[2]);
		A = matrix_product(A,vvvi[2]);
		vvi_a = inverse(vvvi[1])[1];
		vvi_ker = matrix_product(vvi_ker, vvi_a);
		A = matrix_product(vvvi[1], A);
		//show_matrix(mat_add(mat_scalar_prod(vvi_im,-1), matrix_product(vvi_ker,A)));
		//show_matrix(A);
		//show_matrix(vvi_ker);
		cout << "rank_ker == " << rank_ker << endl;
		text << "rank_ker == " << rank_ker << endl;
		cout << "rank_im == " << rank_im << endl;
		text << "rank_im == " << rank_im << endl;
		representative.resize(rank_ker - rank_im);
		for (int i = 0; i < rank_ker - rank_im; i++) {
			const int n = vvi_ker.size();
			const int m = vvi_ker[0].size();
			for (int j = 0; j < n; j++) {
				representative[i].push_back(vvi_ker[j][m - i - 1]);
			}
		}
		for (int i = 0; i < representative.size(); i++) {
			cout << endl << "representative " << i << endl;
			text << endl << "representative " << i << endl;
			for (int j = 0; j < representative[i].size(); j++) {
				if (representative[i][j] != 0) {
					cout << representative[i][j] << "," << j << endl;
					text << representative[i][j] << "," << j << endl;
				}
			}
		}
	}
	if (true) {
		assert(representative.size() != 0);
		vi = {};
		vi.resize(representative[0].size(), 0);
		all_representative.resize(pow(2,representative.size()), vi);
		for (int i = 0; i < pow(2, representative.size()); i++) {
			for (int j = 0; j < representative.size(); j++) {
				if (i & (1 << j)) {
					for (int k = 0; k < all_representative.size(); k++) {
						all_representative[i][k] += representative[j][k];
					}
				}
			}
		}
		vvi_a = ch.base[4];
		vvi_b = ch.base[4];
		for (int k = 0; k < ch.rank[4]; k++) {
			vvi_a[k].erase(vvi_a[k].begin() + 3, vvi_a[k].begin() + 4);
			vvi_b[k].erase(vvi_b[k].begin(), vvi_b[k].begin() + 1);
		}
		for (int i = 0; i < all_representative.size(); i++) {
			cout << "The value of all_representative[" << i << "] is ";
			text << "The value of all_representative[" << i << "] is ";
			sum = 0;
			for (int j = 0; j < ch.rank[4]; j++) {
				for (int k = 0; k < all_representative[i].size(); k++) {
					for (int l = 0; l < all_representative[i].size(); l++) {
						if (vvi_a[j] == ch.base[2][k] && vvi_b[j] == ch.base[2][l]) {
							sum = sum + vg_orinent[j].coefficient * all_representative[i][k] * all_representative[i][l];
 						}
					}
				}
			}
			cout << sum << endl;
			text << sum << endl;
		}
	}
}