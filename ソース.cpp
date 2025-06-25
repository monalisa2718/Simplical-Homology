#include<iostream>
#include<string>
#include<vector>

using namespace std;
class generator_of_module {
public:
	vector<string> name;
	int coefficient=0;
};
class simplical_chain {
public:
	vector<vector<string>> simplical_complex;
	vector<vector<string>> simplical_precomplex;
	vector<int> rank;
	vector<vector<vector<string>>> base;
	vector<vector<int>> boundary_matrix;
	vector<vector<int>> coboundary_matrix;
	void set_rank_base() {
		int max_degree = 0;
		for (int i = 0; i < simplical_complex.size(); i++) {
			if (max_degree < simplical_complex[i].size()) {
				max_degree = simplical_complex[i].size();
			}
		}
		for (int i = 1; i <= max_degree; i++) {
			vector<vector<string>> vvs;
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
				}
				cout << "},";
			}
			cout << endl;
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
		//この関数でn-1次のコチェインからn次のチェインへのコバウンダリー写像を作ります。
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
		boundary_matrix.resize(coboundary_matrix[0].size(),vi);
		for (int i = 0; i < coboundary_matrix.size(); i++) {
			for (int j = 0; j < coboundary_matrix[0].size(); j++) {
				boundary_matrix[j][i] = coboundary_matrix[i][j];
			}
		}
	}
	void complexed() {
		simplical_complex.resize(0);
		for (vector<string> vs : simplical_precomplex) {
			for (int bit = 0; bit < (1 << vs.size()); bit++) {
				vector<string> vs_a;
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
	vector<generator_of_module> boundary(vector<string> element_of_chain) {
		vector<generator_of_module> vg;
		element_of_chain.size();
		for (int i = 0; i < element_of_chain.size(); i++) {
			vg.resize(vg.size() + 1);
			for (int j = 0; j < element_of_chain.size(); j++) {
				if (i != j) {
					vg[i].name.push_back(element_of_chain[j]);
				}
			}
			vg[i].coefficient = pow(-1,i);
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
	vi.resize(matrix_b[0].size(),0);
	vvi.resize(matrix_a.size(),vi);
	for (int i = 0; i < vvi.size(); i++) {
		for (int j = 0; j < vvi[0].size(); j++) {
			for (int k = 0; k < matrix_a[0].size(); k++) {
				vvi[i][j] += matrix_a[i][k] * matrix_b[k][j];
			}
		}
	}
	return vvi;
}
vector<vector<int>> add_column(vector<vector<int>>*matrix, int n, int m, int l) {//n列目をm列目にl倍して加える関数。
		vector<vector<int>> vvi = id((*matrix)[0].size());
		vvi[n][m] += l;
		(*matrix) = matrix_product((*matrix), vvi);
		return vvi;
	}
vector<vector<int>> add_row(vector<vector<int>>*matrix, int n, int m, int l) {//n行目をm行目にl倍して加える関数。
		vector<vector<int>> vvi = id((*matrix).size());
		vvi[m][n] += l;
		(*matrix) = matrix_product(vvi, *matrix);
		return vvi;
	}
vector<vector<int>> exchange_column(vector<vector<int>>*ptmatrix, int n, int m) {
		vector<vector<int>> vvi = id((*ptmatrix)[0].size());
		vvi[n][n] = 0;
		vvi[m][m] = 0;
		vvi[m][n] = 1;
		vvi[n][m] = 1;
		(*ptmatrix) = matrix_product((*ptmatrix), vvi);
		return vvi;
	}
vector<vector<int>> exchange_row(vector<vector<int>>*ptmatrix, int n, int m) {
		vector<vector<int>> vvi = id((*ptmatrix).size());
		vvi[n][n] = 0;
		vvi[m][m] = 0;
		vvi[m][n] = 1;
		vvi[n][m] = 1;
		(*ptmatrix) = matrix_product(vvi, (*ptmatrix));
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
		cout << search_number << endl;
		int find_one[3] = { 0,0,0 };
		for (int i = search_number; i < matrix.size(); i++) {
			find_one[2] = 0;
			for (int j = search_number; j < matrix[0].size(); j++) {
				if (vvvi[0][i][j] == 1) {
					find_one[0] = i;
					find_one[1] = j;
					find_one[2]++;
					break;
				}
			}
			if (find_one[2] != 0)break;
		}
		if (find_one[2] != 0) {
			vvvi[1] = matrix_product(exchange_row(&(vvvi[0]), search_number, find_one[0]), vvvi[1]);
			vvvi[2] = matrix_product(vvvi[2], exchange_column(&(vvvi[0]), search_number, find_one[1]));
			for (int i = search_number + 1; i < vvvi[0].size(); i++) {
				vvvi[1] = matrix_product(add_row(&(vvvi[0]), search_number, i, -1 * vvvi[0][i][search_number]), vvvi[1]);
			}
			for (int i = search_number + 1; i < vvvi[0][0].size(); i++) {
				vvvi[2] = matrix_product(vvvi[2], add_column(&(vvvi[0]), search_number, i, -1 * vvvi[0][search_number][i]));
			}
		}
	}
	
	return vvvi;
}

int main() {
	simplical_chain ch;
	ch.simplical_precomplex = { {"0","1","2","3","4"},{"1","2","3","4","5"},{"2","3","4","5","6"},{"3","4","5","6","7"},{"4","5","6","7","8"},{"5","6","7","8","9"},{"6","7","8","9","a"},{"7","8","9","a","b"},{"0","1","9","a","b"},{"0","1","2","a","b"},{"0","1","2","3","b"}};//頂点の順序を守って表記しないとおかしなことになる。
	ch.complexed();
	ch.set_rank_base();
	ch.show_base();
	ch.set_coboundary_matrix(3);

	vector<vector<int>> vvi_a = { {4,2},{3,1},{1,6} };
	vector<vector<int>> vvi_b = { {1,2,3} };
	vector<vector<int>> vvi_c;
	vector<vector<vector<int>>> vvvi;
	cout << endl << "delta 2" << endl;
	show_matrix(ch.coboundary_matrix);
	vvvi = Smith_normalize(ch.coboundary_matrix);
	show_matrix(vvvi[0]);
	show_matrix(vvvi[1]);
	show_matrix(vvvi[2]);
	cout << endl << "deta 1" << endl;
	ch.set_coboundary_matrix(2);
	show_matrix(ch.coboundary_matrix);
	vvvi.resize(0);
	vvvi = Smith_normalize(ch.coboundary_matrix);
	show_matrix(vvvi[0]);
	show_matrix(vvvi[1]);
	show_matrix(vvvi[2]);

}