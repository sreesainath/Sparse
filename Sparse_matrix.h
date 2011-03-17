#include<iostream>
using namespace std;

class Sparse_matrix
{

private:
//The data types are specified int initially for simplicity. These will be converted to Template format soon.
	int num_dims;               // No. of dimensions or Dimensionality
	int *shape;                 // Shape of the matrix. For e.g. 3 X 2 X 2
	int num_nzvals;             // No. of non-zero values in the sparse matrix
	int *nz_vals;               // Array of non-zero values
	int **ptr;                  // Pointer to the values(compressed format)
	int **indices;  	    // Indices for non- zero values

public:
	Sparse_matrix(int dim, int num_val, int *val, int *shape, int ** indices)
	{
	// Initializations
		num_dims = dim;
		num_nzvals = num_val;
		nz_vals = new int[num_nzvals];
		for(int i = 0; i < num_nzvals; i++)
			nz_vals[i] = val[i];
		shape = new int[num_dims];
		for(int i = 0; i < num_dims; i++)
			shape[i] = shape[i]; 
		ptr = new int*[num_dims];
		for(int i =0; i < num_dims; i++)
			ptr[i] = NULL;
		indices = new int*[num_nzvals];
		for(int i =0; i < num_dims; i++)
			indices[i] = new int[num_dims];
		
		for(int i =0; i < num_nzvals; i++)
			for(int j = 0; j < num_dims; j++)
			indices[i][j] = indices[i][j];
	
	} 
	Sparse_matrix()
	{
	// Initialization of an empty Sparse_matrix object
		num_dims = 0;
		num_nzvals = 0;
		nz_vals = new int[num_nzvals];
		ptr = new int*[num_dims];
		for(int i =0; i < num_dims; i++)
			ptr[i] = NULL;
		indices = new int*[num_nzvals];
		for(int i =0; i < num_dims; i++)
			indices[i] = new int[num_dims];
	}
	~Sparse_matrix()
	{
	//Destructor
		delete[] nz_vals;
	}
	void printdata();
//	Sparse_matrix addition(Sparse_matrix);
	Sparse_matrix square();
//	void square();
//	void sq_root();
	int cumprod();
	int cumsum();
};

void Sparse_matrix :: printdata()
{
// Printing values to moniter
	cout<< "\nDimensionality: " << num_dims;
	cout<< "\nNo. of non-zero Values: " << num_nzvals;
	cout<< "\nValues: ";
	for(int i = 0; i < num_nzvals; i++)
		cout<< nz_vals[i] << "\t";
}

/*Sparse_matrix Sparse_matrix :: square()
{
// Square of the matrix
	Sparse_matrix a(num_dims,num_nzvals,nz_vals,shape,indices);
	for(int i = 0; i < num_nzvals; i++)
		a.nz_vals[i] = nz_vals[i] * nz_vals[i];
	return a;
}
/*void Sparse_matrix :: square()
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = nz_vals[i] * nz_vals[i];
}*/

int Sparse_matrix :: cumprod()
{
// Product of all non-zero elements
	int prod = 1;
	for(int i = 0; i < num_nzvals; i++)
		prod *= nz_vals[i];
	return prod;
}
int Sparse_matrix :: cumsum()
{
// Product of all non-zero elements
	int sum = 0;
	for(int i = 0; i < num_nzvals; i++)
		sum += nz_vals[i];
	return sum;
}
/*Sparse_matrix Sparse_matrix :: addition(Sparse_matrix b)
{
//Matrix addition
	// Step 1: Compatibility Check
	if(num_dims != b.num_dims) 
	{
		cout<<"\n Matrices Incompatible for addition";
		return;
	}
	else
	{ // This block will be removed once != operator is overloaded
		for(int i = 0; i < num_dims; i++)
		{
			if(shape[i] != b.shape[i])
			{
				cout<<"\nMatrices Incompatible for addition";
				return;
			}
		}
	}
}*/

