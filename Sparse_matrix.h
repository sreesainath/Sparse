#include<iostream>
#include<math.h>
using namespace std;

class Sparse_matrix
{

private:
//The data types are specified int initially for simplicity. These will be converted to Template format soon. Values in the sparse matrix are considered to be real. Assuming for now that the matricces are index sorted.

	int num_dims;               // No. of dimensions or Dimensionality
	int *shape;                 // Shape of the matrix. For e.g. 3 X 2 X 2
	int num_nzvals;             // No. of non-zero values in the sparse matrix
	double *nz_vals;            // Array of non-zero values
	int **ptr;                  // Pointer to the values(compressed format)
	int **indices;  	    // Indices for non- zero values

public:
// Member functions
	void printdata();
	// Arithmetic Operations
	double cumprod();	     // To be changed as a single function that may take axis as arguement
	double cumsum();	     // To be changed as a single function that may take axis as arguement
	double max();
	double min();
	double mean();
	Sparse_matrix square();
	Sparse_matrix sq_root();
//	void square();
//	Sparse_matrix addition(Sparse_matrix);

Sparse_matrix(int dim, int num_val, double *val, int *shape, int ** indices)
{
// Initializations in case values are known. Otherwise default constructor is invoked.
	num_dims = dim;
	num_nzvals = num_val;
	nz_vals = new double[num_nzvals];
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = val[i];
	shape = new int[num_dims];
	for(int i = 0; i < num_dims; i++)
		shape[i] = shape[i]; 
	ptr = new int*[num_dims];
	for(int i =0; i < num_dims; i++)
		ptr[i] = NULL;
	indices = new int*[num_nzvals];
	for(int i =0; i < num_nzvals; i++)
		indices[i] = new int[num_dims];
	for(int i =0; i < num_nzvals; i++)
		for(int j = 0; j < num_dims; j++)
		indices[i][j] = indices[i][j];
} 

~Sparse_matrix()
{
//Destructor - tbc
	delete[] nz_vals;
}

};


void Sparse_matrix :: printdata()		// Printing values to moniter
{
	cout<< "\nDimensionality: " << num_dims;
	cout<< "\nNo. of non-zero Values: " << num_nzvals;
	cout<< "\nValues: ";
	for(int i = 0; i < num_nzvals; i++)
		cout<< nz_vals[i] << "\t";
}

Sparse_matrix Sparse_matrix :: square()		// Square of the matrix
{
	Sparse_matrix a(num_dims,num_nzvals,nz_vals,shape,indices);
	for(int i = 0; i < num_nzvals; i++)
	{
		a.nz_vals[i] = nz_vals[i] * nz_vals[i];
		cout<<"\n" << a.nz_vals[i];
	}
	return a;
}

Sparse_matrix Sparse_matrix :: sq_root()		// Square root of the matrix
{
	Sparse_matrix a(num_dims,num_nzvals,nz_vals,shape,indices);
	for(int i = 0; i < num_nzvals; i++)
		a.nz_vals[i] = sqrt(nz_vals[i]);
	return a;
}

double Sparse_matrix :: cumprod()			// Product of all non-zero elements - Cumulative Product
{
	double prod = 1;
	for(int i = 0; i < num_nzvals; i++)
		prod *= nz_vals[i];
	return prod;
}

double Sparse_matrix :: cumsum()	
{		
	double sum = 0;
	for(int i = 0; i < num_nzvals; i++)
		sum += nz_vals[i];
	return sum;
}

double Sparse_matrix :: max()		// Maximum value in the matrix , will be extended to cover axes in the next revision
{
	double max = nz_vals[0];
	for(int i = 1; i < num_nzvals; i++)
	{
		if(nz_vals[i] > max)
			max = nz_vals[i];
	}
	return max;
}

double Sparse_matrix :: min()		// Miniimum value in the matrix , will be extended to cover axes in the next revision
{
	double min = nz_vals[0];
	for(int i = 1; i < num_nzvals; i++)
	{
		if(nz_vals[i] < min)
			min = nz_vals[i];
	}
	return min;
}

double Sparse_matrix :: mean()		// Mean value in the matrix , will be extended to cover axes in the next revision
{
	double mean = 0.0;
	for(int i = 0; i < num_nzvals; i++)
	{
		mean += nz_vals[i];
	}
	mean = mean / num_nzvals;
	return mean;
}

/*void Sparse_matrix :: square()		// Same matrix is used to store the squared value
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = nz_vals[i] * nz_vals[i];
}*/


Sparse_matrix Sparse_matrix :: addition(Sparse_matrix b)	//Matrix addition - index-wise
{
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
		Sparse_matrix a;
		int nz = 0;
		if(num_nzvals > b.num_nzvals)
		{
			for(int i = 0; i < num_nzvals ; i++)
				for(int k = 0; k < num_dims; k++)
				{
					if(indices[i][k] == b.indices[i][k])
						continue;
					else
						nz++;
				}
		}
		a.num_nzvals = num_nzvals + nz;
	}
}

