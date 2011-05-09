#include<iostream>
#include<math.h>
using namespace std;

class Sparse_matrix
{

private:
//The data types are specified int initially for simplicity. These will be converted to Template format soon. Values in the sparse matrix are considered to be real. Assuming for now that the matrices are index sorted.
// Function names have been mirrored as much as possible to those in numpy.ndarray

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
	double max();		     // Maximum value
	int* argmax();		     // Indices of the maximum value
	double min();                // Minimum value
	int* argmin();               // Indices of the minimum value
	double mean();               // Mean value
	void sort();                 // Sort the matrix
	double var();		     // Variance
	double std();                // Standard Deviation
//	int** argsort();
	void __neg__();              // Negative 
	void __pos__();              // Positive
	void __abs__();              // Absolte
	double* ptp(int axis = 0);   // Maxima and minima along an axis
	void clip(double min, double max);        // Clip the values to an interval
	void swapaxes(int axis1, int axis2);      // swap values in two axes
	void flatten();                           // Flatten a matrix
	double getitem(int *index);               // Get an item by index
	void setitem(int *index, double value);   // Set the value of an item
//	Sparse_matrix square(); 
	Sparse_matrix sq_root();                  // Square root
	void square();                            // Square
//	Sparse_matrix addition(Sparse_matrix);

Sparse_matrix(int dim, int num_val, double *val, int *shape, int **index)
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
	ptr = new int* [num_nzvals];
	for(int i =0; i < num_nzvals; i++)
		ptr[i] = new int[num_dims];
	indices = new int*[num_nzvals];
	for(int i =0; i < num_nzvals; i++)
		indices[i] = new int[num_dims];
	for(int i =0; i < num_nzvals; i++)
		for(int j = 0; j < num_dims; j++)
		indices[i][j] = index[i][j];
	for(int i = 0; i < num_nzvals; i++)
		ptr[i][0] = indices[i][1];
	for(int j =1; j < num_dims ; j++)
		ptr[0][j] = 0;
	int j = 1,l = 1,diff = 0;
	for(int i = 1; i < num_nzvals; i++)
		{
			diff = indices[i][0] - indices[i-1][0];
			if(diff == 0)
				continue;
			else
			{
				ptr[l][j] = i;
				l++;
			}
		}
	j = 2;
	while(j < num_dims)
	{
		l = 1;
		for(int i = 1; i < num_nzvals; i++)
                {
                        diff = indices[i][j] - indices[i-1][j];
                        if(diff == 0)
                                continue;
                        else
                        {
                                ptr[l][j] = i;
                                l++;
                        }
                }
		j++;
        }
} 

~Sparse_matrix()
{
//Destructor - tbc
	delete[] nz_vals;
}

};

// Member Functions
//1. printdata()
void Sparse_matrix :: printdata()		// Printing values to moniter
{
	cout<< "\nDimensionality: " << num_dims;
	cout<< "\nNo. of non-zero Values: " << num_nzvals;
	cout<< "\nValues: ";
	for(int i = 0; i < num_nzvals; i++)
		cout<< nz_vals[i] << "\t";
	cout<<"\nIndices:";
	for(int i = 0; i < num_nzvals; i++)
	{
		cout<<"\n";
		for(int j = 0; j < num_dims; j++)
			cout<<"\t"<< indices[i][j];
	}
	cout<<"\n Compressed format:";
	for(int i = 0; i < num_nzvals; i++)
        {
                cout<<"\n";
                for(int j = 0; j < num_dims; j++)
                        cout<<"\t"<< ptr[i][j];
        }

}

/*
//2. square()
Sparse_matrix Sparse_matrix :: square()		// Square of the matrix - the same matrix stores the output values
{
	Sparse_matrix a(num_dims,num_nzvals,nz_vals,shape,indices);
	for(int i = 0; i < num_nzvals; i++)
	{
		a.nz_vals[i] = nz_vals[i] * nz_vals[i];
		cout<<"\n" << a.nz_vals[i];
	}
	return a;
}
*/

//3. sq_root()
Sparse_matrix Sparse_matrix :: sq_root()		// Square root of the matrix
{
	Sparse_matrix a(num_dims,num_nzvals,nz_vals,shape,indices);
	for(int i = 0; i < num_nzvals; i++)
		a.nz_vals[i] = sqrt(nz_vals[i]);
	return a;
}

//4. Negative
void Sparse_matrix :: __neg__()           // -A
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = -1 * nz_vals[i];
}

//5. Positive
void Sparse_matrix :: __pos__()             // +A
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = +1 * nz_vals[i];
}

//6. Absolute value
void Sparse_matrix :: __abs__()                  // Absolute value of each element in A
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = fabs(nz_vals[i]);
}
	
//7. Cumulative product
double Sparse_matrix :: cumprod()			// Product of all non-zero elements - Cumulative Product
{
	double prod = 1;
	for(int i = 0; i < num_nzvals; i++)
		prod *= nz_vals[i];
	return prod;
}

//8. Cumulative sum
double Sparse_matrix :: cumsum()	// Cumulative Sum
{		
	double sum = 0;
	for(int i = 0; i < num_nzvals; i++)
		sum += nz_vals[i];
	return sum;
}

//9. Peak to peak along an axis
double* Sparse_matrix :: ptp(int axis = 0)
{
	double minima = min(axis);
	double maxima = max(axis);
	double peaks[2] = {maxima,minima};
	return peaks;
}

//10. Product along an axis
/*
Sparse_matrix Sparse_matrix :: prod(int axis = 0)
{
	double **shape1;
	
	for(int i = 1; i < num_nzvals[i]; i++)
	{
		if(indices[i][0] == indices[i-1][0])	
			prod[i] *= nz_vals[i];
		else
			prod[i] = nz_vals[i];
	}
		
}
*/

//11. Maximum value along an axis
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

//12. Argmax
int* Sparse_matrix :: argmax()          // Returns the indices of the maximum value
{
	int *argsmax;
	argsmax = new int[num_dims];
	argsmax = indices[0];
	double max = nz_vals[0];
        for(int i = 1; i < num_nzvals; i++)
        {
                if(nz_vals[i] > max)
		{
                        max = nz_vals[i];
			argsmax = indices[i];
		}
        }
        return argsmax;
}

//13. Minimum value
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

//14. Argmin
int* Sparse_matrix :: argmin()         // Returns the indices of the minimum value
{
        int *argsmin;
        argsmin = new int[num_dims];
        argsmin = indices[0];
        double min = nz_vals[0];
        for(int i = 1; i < num_nzvals; i++)
        	if(nz_vals[i] < min)
                	argsmin = indices[i];
       	return argsmin;
}

//15. Mean value
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

double Sparse_matrix :: var()
{

	double mean = mean();
	double var[num_nzvals], variance = 0.0;
	for(int i = 0; i < num_nzvals; i++)
	{
		var[i] = mean - nz_vals[i];
		variance = variance + (var[i]*var[i]);
	}
	variance = variance / num_nzvals;
	return variance;
}

double Sparse_matrix :: std()
{
	double var = var();
	double std = sqrt(var);
	return std;
}

//16. Square
void Sparse_matrix :: square()		// Same matrix is used to store the squared value
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = nz_vals[i] * nz_vals[i];
}

//17. Sum of two matrices
/*
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
*/

//18. Sort the matrix
void Sparse_matrix :: sort()
{
	double temp = 0.0000;
	for(int i = 0; i < num_nzvals; i++)  // Sorting only changes the values in the matrix, not the positions. Only ascending considered for now.
	for(int j = 0; j < num_nzvals; j++)
	{
		if(nz_vals[i] < nz_vals[j])
		{	
			temp = nz_vals[i];
			nz_vals[i] = nz_vals[j];
			nz_vals[j] = temp;
		}
			
	}
}

//19. Get an item
double Sparse_matrix :: getitem(int *index)  // Get the item at a specific location
{
	int flag = 0;
	for(int i = 0; i < num_nzvals; i++)
	{
		for(int j = 0; j < num_dims; j++)
		{
			if(indices[i][j] == index[j])
			{
				flag = 1;
				continue;
			}
			else
			{
				flag = 0;
				break;
			}
		}
		if(flag == 1)
			return nz_vals[i];
		else
			continue;
	}
	return 0;
}


//20. Set value for an item
void Sparse_matrix :: setitem(int *index, double value)   // Set the value at a specific location
{
	int flag = 0;
        for(int i = 0; i < num_nzvals; i++)
        {
                for(int j = 0; j < num_dims; j++)
                {
                        if(indices[i][j] == index[j])
                        {
                                flag = 1;
                                continue;
                        }
                        else
                        {
                                flag = 0;
                                break;
                        }
                }
                if(flag == 1)
                        nz_vals[i] = value;
                else
                        continue;
        }
	if(flag == 0)
	{
		num_nzvals += 1;
		nz_vals[i] = value;
		for(int j = 0; j <  num_dims; j++)
			indices[i][j] = index[j];
	}

}

//21. Swap axes
void Sparse_matrix :: swapaxes(int axis1 = 0, int axis2 = 1)    // Swap two axes of the matrix. In case of a 2D matrix, transpose of the matrix
{
	int temp = shape[axis1];
	shape[axis1] = shape[axis2];
	shape[axis2] = temp;
	for(int i = 0; i < num_nzvals; i++)
	{
		int temp = indices[i][axis1];
		indices[i][axis1] = indices[i][axis2];
		indices[i][axis2] = temp;
	}
	for(int i = 0; i < num_nzvals; i++)
                ptr[i][0] = indices[i][1];
        for(int j =1; j < num_dims ; j++)
                ptr[0][j] = 0;
        int j = 1,l = 1,diff = 0;
        for(int i = 1; i < num_nzvals; i++)
                {
                        diff = indices[i][0] - indices[i-1][0];
                        if(diff == 0)
                                continue;
                        else
                        {
                                ptr[l][j] = i;
                                l++;
                        }
                }
        j = 2;
        while(j < num_dims)
        {
                l = 1;
                for(int i = 1; i < num_nzvals; i++)
                {
                        diff = indices[i][j] - indices[i-1][j];
                        if(diff == 0)
                                continue;
                        else
                        {
                                ptr[l][j] = i;
                                l++;
                        }
                }
                j++;
        }
}

//22. clip 
void Sparse_matrix :: clip(double min, double max)    // Clip the matrix values to a value range
{
	for(int i = 0; i < num_nzvals; i++)
	{
		if(nz_vals[i] < min)
			nz_vals[i] = min;
		else
		{
			if(nz_vals[i] > max)
				nz_vals[i] = max;
		}
	}
}

//23. Flatten
void Sparse_matrix :: flatten()                      // Convert the matrix into a 1D array
{
	int vals = 1;
	for(int i = 0; i < num_dims; i++)
		vals *= shape[i];
	shape[0] = 1; shape[2] = vals;
	int index[num_nzvals][2];
	for(int i = 0; i < num_nzvals; i++)
	{
		if(indices[i+1][0] - indices[i][0] >= 0)
		{
			index[i][0] = 0;
			index[i][1] = 6 * indices[i][0];
			for(int j = 1; j < num_dims; j++)
				index[i][1] = index[i][1] + indices[i][j];
		}
		else
		{
			index[i][0] = 0;
			index[i][1] = index[i-1][1] + (6*indices[i][0])+ indices[i][1];
		}
		indices[i][0] = 0;		
		indices[i][1] = index[i][1];
	}
	
}
/*
Sparse_matrix Sparse_matrix :: take(int **index, axis = 0)
{
	
}
*/
/*
Sparse_matrix Sparse_matris :: put(int **index, double *values)
{
}
*/
/*
void fill(double x)
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = x;
}
*/

