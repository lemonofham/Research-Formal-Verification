#include <algorithm>
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>
#include <tuple>
#include <cmath>

using namespace std;


typedef tuple<int, vector<int>> int_vector_tuple;
void block_indice_values(int i, int block_length, int grid_length, vector<int>& all_block_indices);
bool clause(int i, int grid_length, int block_length, const vector<int>& literals, const vector<int>& grid, const vector<int>& all_block_indices, bool extra=false, const vector<int>& grid_for_clause={});
void update_literals(int grid_length, int block_length, vector<int>& literals, int update_variable, vector<int>& grid, vector<int>& number_of_false, const vector<int>& all_block_indices);
void set_literals(int i, int grid_length, int block_length, vector<int>& literals, vector<int>& grid, vector<int>& number_of_false, const vector<int>& all_block_indices);
int_vector_tuple depth_first_search(int grid_length, int block_length, vector<int>& grid, vector<int>& literals, vector<int>& number_of_false, const vector<int>& all_block_indices, bool extra=false, const vector<int>& grid_for_clause={});
int_vector_tuple sudoku_sat_solver(vector<int>& grid);
void initial_setup(int grid_length, vector<int>& grid, const vector<int>& all_block_indices);
void print_sudoku(int grid_length, const vector<int>& grid, int start_index=0);

int main()
{
	int n = 9;
	// int n;
	// cin>>n;
	int grid_size = n * n;
	vector<int> grid(grid_size, 0);
	grid = {0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 2, 0, 3, 0, 0, 0, 0, 4, 1, 7, 0, 0, 6, 0, 9, 5, 0, 3, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 7, 0, 6, 1, 0, 0, 7, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 5, 0, 0, 0, 9, 7, 0, 0, 0, 0, 9, 3, 8, 0};
	// for(int i = 0; i < grid_size; i++)
	// 	cin>>grid[i];
	auto start = std::chrono::high_resolution_clock::now();
	int final_result;
	vector<int> returned_grid;
	tie(final_result, returned_grid) = sudoku_sat_solver(grid);
	if(final_result == 0)
		cout<<"Wrong Sudoku!!"<<endl;
	else if(final_result == 1)
	{
		cout<<"Multiple Solutions Exist for this Sudoku!!\nTwo such are:"<<endl;
		print_sudoku(n, returned_grid);
		cout<<endl;
		print_sudoku(n, returned_grid, grid_size);
	}
	else
	{
		cout<<"Solution Exists!!"<<endl;
		print_sudoku(n, returned_grid);
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto d = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
	if(d < 1000)
		cout<<fixed<<setprecision(2)<<"Time taken to run: "<<d<<"μs!"<<endl;
	else if(d > 1000 && d < 1000000)
		cout<<fixed<<setprecision(2)<<"Time taken to run: "<<((double) d / 1000)<<"ms!"<<endl;
	else if(d > 1000000)
	{
		double d_seconds = (double) d / 1000000;
		if(d_seconds < 60)
			cout<<fixed<<setprecision(2)<<"Time taken to run: "<<(d_seconds)<<"s!"<<endl;
		if(d_seconds > 60 && d_seconds < 3600)
			cout<<fixed<<setprecision(2)<<"Time taken to run: "<<((double) d_seconds / 60)<<"min(s)!"<<endl;
		if(d_seconds < 60)
			cout<<fixed<<setprecision(2)<<"Time taken to run: "<<((double) d_seconds / 3600)<<"hr(s)!"<<endl;
	}
}

int_vector_tuple sudoku_sat_solver(vector<int>& grid)
{
	int satisfiable = 0;
	int once_satisfiable = 0;
	int twice_satisfiable;
	int total_cells = grid.size();
	int grid_length = (int) sqrt(total_cells);
	int block_length = (int) sqrt(grid_length);
	vector<int> satisfied_grid;
	vector<int> another_satisfied_grid;
	vector<int> returning_grid(2 * total_cells, 0);
	vector<int> literals(grid_length * total_cells, 0);
	vector<int> number_of_false(total_cells, 0);
	vector<int> temp_grid(grid);
	vector<int> all_block_indices(total_cells, 0);
	for(int i = 0; i < grid_length; i++)
		block_indice_values((block_length * ((i % block_length) + grid_length * (i / block_length))), block_length, grid_length, all_block_indices);
	// initial_setup(grid_length, temp_grid, all_block_indices);
	for(int i = 0; i < total_cells; i++)
		if(temp_grid[i] != 0)
			set_literals(i, grid_length, block_length, literals, temp_grid, number_of_false, all_block_indices);	
	vector<int> temp_temp_grid(temp_grid);
	vector<int> temp_literals(literals);
	vector<int> temp_number_of_false(number_of_false);
	tie(satisfiable, satisfied_grid) = depth_first_search(grid_length, block_length, temp_grid, literals, number_of_false, all_block_indices);
	if(satisfiable)
	{
		copy(satisfied_grid.begin(), satisfied_grid.end(), returning_grid.begin());
		temp_grid = temp_temp_grid;
		literals = temp_literals;
		number_of_false = temp_number_of_false;
		tie(twice_satisfiable, another_satisfied_grid) = depth_first_search(grid_length, block_length, temp_grid, literals, number_of_false, all_block_indices, true, satisfied_grid);
		once_satisfiable = !twice_satisfiable;
		if(twice_satisfiable)
			copy(another_satisfied_grid.begin(), another_satisfied_grid.end(), returning_grid.begin() + total_cells);
	}
	return make_tuple(satisfiable + once_satisfiable, returning_grid);
}

int_vector_tuple depth_first_search(int grid_length, int block_length, vector<int>& grid, vector<int>& literals, vector<int>& number_of_false, const vector<int>& all_block_indices, bool extra, const vector<int>& grid_for_clause)
{
	bool solved = true;
	for(int i = 0; i < grid.size(); i++)
		solved = solved & clause(i, grid_length, block_length, literals, grid, all_block_indices, extra, grid_for_clause);
	if(solved)
		return make_tuple(solved, grid); 
	vector<int> temp_grid(grid);
	vector<int> temp_literals(literals);
	vector<int> returning_grid(grid.size(), 0);
	for(int i = 0; i < grid_length * grid.size(); i++)
	{
		if(literals[i] == 0)
		{
			temp_grid[i / grid_length] = (i % grid_length) + 1;
			set_literals(i / grid_length, grid_length, block_length, temp_literals, temp_grid, number_of_false, all_block_indices);
			tie(solved, returning_grid) = depth_first_search(grid_length, block_length, temp_grid, temp_literals, number_of_false, all_block_indices, extra, grid_for_clause);
			if(solved)
				return make_tuple(solved, returning_grid);
			else
			{
				temp_grid = grid;
				temp_literals = literals;
				temp_literals[i] = 2;
				tie(solved, returning_grid) = depth_first_search(grid_length, block_length, temp_grid, temp_literals, number_of_false, all_block_indices, extra, grid_for_clause);
				return make_tuple(solved, returning_grid);
			}
			break;
		}
	}
	return make_tuple(false, returning_grid);
}

void set_literals(int i, int grid_length, int block_length, vector<int>& literals, vector<int>& grid, vector<int>& number_of_false, const vector<int>& all_block_indices)
{
	for(int j = 0; j < grid_length; j++)
	{
		if(grid[i] - 1 == j)
		{
			literals[grid_length * i + j] = 1;
			continue;
		}
		literals[grid_length * i + j] = 2;
		number_of_false[i] += 1;
	}
	int row_number = i / grid_length;
	int column_number = i % grid_length;
	int start = grid_length * (column_number / block_length + block_length * (row_number / block_length));
	for(int j = 0; j < grid_length; j++)
	{
		if(j == column_number || grid[(grid_length * row_number + j)] != 0)
			continue;
		int update_variable = (grid_length * (grid_length * row_number + j)) + grid[i] - 1;
		int updating = update_variable / grid_length;
		if(literals[update_variable] != 0)
			continue;
		literals[update_variable] = 2;
		number_of_false[updating] += 1;
		if(number_of_false[updating] == grid_length - 1)
			update_literals(grid_length, block_length, literals, updating, grid, number_of_false, all_block_indices);
	}
	for(int j = 0; j < grid_length; j++)
	{
		if(j == row_number || grid[(grid_length * j + column_number)] != 0)
			continue;
		int update_variable = (grid_length * (grid_length * j + column_number)) + grid[i] - 1;
		int updating = update_variable / grid_length;
		if(literals[update_variable] != 0)
			continue;
		literals[update_variable] = 2;
		number_of_false[updating] += 1;
		if(number_of_false[updating] == grid_length - 1)
			update_literals(grid_length, block_length, literals, updating, grid, number_of_false, all_block_indices);
	}
	for(int j = 0; j < grid_length; j++)
	{
		if(all_block_indices[start + j] == i || grid[all_block_indices[start + j]] != 0)
			continue;
		int update_variable = (grid_length * all_block_indices[start + j]) + grid[i] - 1;
		int updating = update_variable / grid_length;
		if(literals[update_variable] != 0)
			continue;
		literals[update_variable] = 2;
		number_of_false[updating] += 1;
		if(number_of_false[updating] == grid_length - 1)
			update_literals(grid_length, block_length, literals, updating, grid, number_of_false, all_block_indices);
	}
	return;
}

void update_literals(int grid_length, int block_length, vector<int>& literals, int updating, vector<int>& grid, vector<int>& number_of_false, const vector<int>& all_block_indices)
{
	int which_true = 0;
	int start = grid_length * updating;
	for(int i = 0; i < grid_length; i++)
		if(literals[start + i] != 2)
		{
			which_true = i + 1;
			break;
		}
	grid[updating] = which_true;
	set_literals(updating, grid_length, block_length, literals, grid, number_of_false, all_block_indices);
	return;
}

bool clause(int i, int grid_length, int block_length, const vector<int>& literals, const vector<int>& grid, const vector<int>& all_block_indices, bool extra, const vector<int>& grid_for_clause)
{
	int total_cells = grid.size();
	int row_number = i / grid_length;
	int column_number = i % grid_length;
	int start = grid_length * (column_number / block_length + block_length * (row_number / block_length));
	bool cell_clause = false;
	bool row_clause = true;
	bool column_clause = true;
	bool block_clause = true;
	bool extra_clause = false;
	if(extra)
	{
		for(int i = 0; i < grid_for_clause.size(); i++)
			extra_clause = extra_clause | !((bool) (literals[grid_length * i + grid_for_clause[i] - 1] % 2));
		if(!extra_clause)
			return false;
	}
	for(int j = 0; j < grid_length; j++)
		cell_clause = cell_clause | ((bool) (literals[grid_length * i + j] % 2));
	if(!cell_clause)
		return false;
	for(int j = 0; j < grid_length; j++)
		for(int k = j + 1; k < grid_length; k++)
			cell_clause = cell_clause & (!((bool) (literals[grid_length * i + j] % 2)) | !((bool) (literals[grid_length * i + k] % 2)));
	if(!cell_clause)
		return false;
	for(int j = 0; j < grid_length; j++)
		for(int k = j + 1; k < grid_length; k++)
			for(int l = 0; l < grid_length; l++)
				row_clause = row_clause & (!((bool) (literals[grid_length * (grid_length * row_number + j) + l] % 2)) | !((bool) (literals[grid_length * (grid_length * row_number + k) + l] % 2)));
	if(!row_clause)
		return false;
	for(int j = 0; j < grid_length; j++)
		for(int k = j + 1; k < grid_length; k++)
			for(int l = 0; l < grid_length; l++)
				column_clause = column_clause & (!((bool) (literals[grid_length * (grid_length * j + column_number) + l] % 2)) | !((bool) (literals[grid_length * (grid_length * k + column_number) + l] % 2)));
	if(!column_clause)
		return false;
	for(int j = 0; j < grid_length; j++)
		for(int k = j + 1; k < grid_length; k++)
			for(int l = 0; l < grid_length; l++)
				block_clause = block_clause & (!((bool) (literals[grid_length * all_block_indices[start + j] + l] % 2)) | !((bool) (literals[grid_length * all_block_indices[start + k] + l] % 2)));
	if(!block_clause)
		return false;
	return true;
}

void block_indice_values(int i, int block_length, int grid_length, vector<int>& all_block_indices)
{
	int row_number = i / grid_length;
	int column_number = i % grid_length;
	for(int j = 0; j < grid_length; j++)
	{
		all_block_indices[i + (grid_length - block_length) * (column_number / block_length) + j] = (j % block_length + grid_length * (j / block_length)) + i - ((column_number % block_length) + grid_length * (row_number % block_length));
	}
	return;
}

void initial_setup(int grid_length, vector<int>& grid, const vector<int>& all_block_indices)
{
	for(int i = 0; i < grid_length; i++)
	{
		vector<int> temp(grid_length, 0);
		int number_of_numbers = 0;
		int free_cell = 0;
		for(int j = 0; j < grid_length; j++)
		{
			if(grid[i * grid_length + j] != 0)
			{
				number_of_numbers += 1;
				temp[grid[i * grid_length + j] - 1] = 1;
			}
			else
				free_cell = i * grid_length + j;
		}
		if(number_of_numbers == grid_length - 1)
		{
			for(int j = 0; j < grid_length; j++)
			{
				if(temp[j] == 0)
				{
					grid[free_cell] = j + 1;
					break;
				}
			}
			initial_setup(grid_length, grid, all_block_indices);
			return;
		}
	}
	for(int i = 0; i < grid_length; i++)
	{
		vector<int> temp(grid_length, 0);
		int number_of_numbers = 0;
		int free_cell = 0;
		for(int j = 0; j < grid_length; j++)
		{
			if(grid[i * grid_length + j] != 0)
			{
				number_of_numbers += 1;
				temp[grid[j * grid_length + i]] = 1;
			}
			else
				free_cell = j * grid_length + i;
		}
		if(number_of_numbers == grid_length - 1)
		{
			for(int j = 0; j < grid_length; j++)
			{
				if(temp[j] == 0)
				{
					grid[free_cell] = j + 1;
					break;
				}
			}
			initial_setup(grid_length, grid, all_block_indices);
			return;
		}
	}
	for(int i = 0; i < grid_length; i++)
	{
		vector<int> temp(grid_length, 0);
		int number_of_numbers = 0;
		int free_cell = 0;
		for(int j = 0; j < grid_length; j++)
		{
			int indice_to_check = all_block_indices[grid_length * i + j];
			if(grid[indice_to_check] != 0)
			{
				number_of_numbers += 1;
				temp[grid[indice_to_check]] = 1;
			}
			else
				free_cell = indice_to_check;
		}
		if(number_of_numbers == grid_length - 1)
		{
			for(int j = 0; j < grid_length; j++)
			{
				if(temp[j] == 0)
				{
					grid[free_cell] = j + 1;
					break;
				}
			}
			initial_setup(grid_length, grid, all_block_indices);
			return;
		}
	}
}

void print_sudoku(int grid_length, const vector<int>& grid, int start_index)
{
	for(int i = 0; i < grid_length; i++)
	{
		for(int j = 0; j < grid_length; j++)
		{
			cout<<grid[grid_length * i + j + start_index];
			if(j != grid_length - 1)
				cout<<" | ";
		}
		cout<<endl;
		if(i != grid_length - 1)
		{
			for(int j = 0; j < grid_length; j++)
				if(j != 0)
					cout<<"--- ";
				else
					cout<<"-- ";
			cout<<endl;
		}
	}
}

// #include <pybind11/pybind11.h>
// namespace py = pybind11;

// int add(int a, int b) {
//     return a + b;
// }

// PYBIND11_MODULE(mycpp, m) {
//     m.def("add", &add, "A function that adds two numbers");
// }