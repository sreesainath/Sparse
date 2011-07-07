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
	double* cumprod();	     // To be changed as a single function that may take axis as arguement
	double* cumsum();	     // To be changed as a single function that may take axis as arguement
	double max();
	int* argmax();
	double min();
	int* argmin();
	double mean();
	void sort();
//	int** argsort();
	void __neg__();
	void __pos__(); 
	void __abs__();
	void clip(double min, double max);
	void swapaxes(int axis1, int axis2);
	void flatten();
	double getitem(int *index);
	void setitem(int *index, double value);
//	Sparse_matrix square();
	Sparse_matrix sq_root();
	void square();
//	Sparse_matrix addition(Sparse_matrix);
	double* dotprod(Sparse_matrix);

	void getshape();
	void addval(double value,int *index);
Sparse_matrix(int dim, int num_val, double *val, int *shapes, int **index)
{
// Initializations in case values are known. Otherwise default constructor is invoked.
	num_dims = dim;
	num_nzvals = num_val;
	nz_vals = new double[num_nzvals];
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = val[i];
	shape = new int[num_dims];
	for(int i = 0; i < num_dims; i++)
		shape[i] = shapes[i]; 
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

void Sparse_matrix :: getshape()
{
	for(int i = 0; i < num_dims; i++)
		cout<<"\t"<< shape[i];
}

/*void Sparse_matrix :: addval(double value, int *index)
{
	++num_nzvals;	
	//nz_vals[num_nzvals] = new double[1];
	nz_vals[num_nzvals] = value;
	indices[num_nzvals] = new int[num_dims];
	for(int j = 0; j < num_dims ; j++)
		indices[num_nzvals][j] = index[j];
	
}*/
void Sparse_matrix :: printdata()		// Printing values to moniter
{
	cout<< "\nDimensionality: " << num_dims;
	cout<<"\nShape of the matrix:";
	for(int i =0; i < num_dims; i++)
		cout<<"\t"<<shape[i];
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
*/

Sparse_matrix Sparse_matrix :: sq_root()		// Square root of the matrix
{
	Sparse_matrix a(num_dims,num_nzvals,nz_vals,shape,indices);
	for(int i = 0; i < num_nzvals; i++)
		a.nz_vals[i] = sqrt(nz_vals[i]);
	return a;
}

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

void Sparse_matrix :: flatten()                      // Convert the matrix into a 1D array
{
	int vals = 1;
	for(int i = 0; i < num_dims; i++)
		vals *= shape[i];
	num_dims = 2;
	shape[0] = 1; shape[2] = vals;
	int index[num_nzvals][2];
	for(int i = 0; i < num_nzvals; i++)
	{
		index[i][0] = 0;
		index[i][1] = indices[i][0]+indices[i][1]+(5 * indices[i][0]);
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
void Sparse_matrix :: nonzero()
{
	for(int i =0; i < num_nzvals; i++)
	{
		cout<<"\n";
                for(int j = 0; j < num_dims; j++)
			cout<<"\t"<<indices[i][j];
	}
}

void fill(double x)
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = x;
}
*/
void Sparse_matrix :: __neg__()           // -A
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = -1 * nz_vals[i];
}

void Sparse_matrix :: __pos__()             // +A
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = +1 * nz_vals[i];
}

void Sparse_matrix :: __abs__()                  // Absolute value of each element in A
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = fabs(nz_vals[i]);
}
	
double*  Sparse_matrix :: cumprod()			// Product of all non-zero elements - Cumulative Product
{
	double * prod;
	prod = new double[num_nzvals];
	prod[0] = nz_vals[0];;
	for(int i = 0; i < num_nzvals; i++)
		prod[i] = prod[i-1] * nz_vals[i];
	return prod;
}

void Sparse_matrix :: swapaxes(int axis1 = 0, int axis2 = 1)    // Swap two axes of the matrix. In case of a matrix, transpose of the matrix
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

/*
Sparse_matrix Sparse_matrix :: prod(int axis = 0)
{
	double prod[i] = 1;
	for(int i = 1; i < num_nzvals[i]; i++)
	{
		if(indices[i][0] == indices[i-1][0])	
			prod[i] *= nz_vals[i];
		else
			prod[i] = nz_vals[i];
	}
		
}
*/
double* Sparse_matrix :: cumsum()	
{		
	double *sum;
	sum = new double[num_nzvals];
	sum[0] = nz_vals[0];;
	for(int i = 0; i < num_nzvals; i++)
		sum[i] = sum[i-1] + nz_vals[i];
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

void Sparse_matrix :: square()		// Same matrix is used to store the squared value
{
	for(int i = 0; i < num_nzvals; i++)
		nz_vals[i] = nz_vals[i] * nz_vals[i];
}

void Sparse_matrix :: sort()
{
	double temp = 0.0000;
	for(int i = 0; i < num_nzvals; i++)  // Sorting only changes the values in the matrix, ot the positions. Only ascending considered for now.
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

double Sparse_matrix :: getitem(int *index)
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

void Sparse_matrix :: setitem(int *index, double value)
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
//	nzvals[i] = value;
//	indices[i] = index;
}
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
			a.num_nzvals = num_nzvals;
		else
			a.um_nzvals = b.num_nzvals;
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
double* Sparse_matrix :: dotprod(Sparse_matrix a) // INCOMPLETE
{
	int endval;
	if(num_nzvals > a.num_nzvals)
		endval = num_nzvals;
	else
		endval = a.num_nzvals;

	double *dots;
	dots = new double[endval];
	for(int i = 0; i < endval; i++)
		dots[i] = 0;

	if(num_dims != a.num_dims)
	{
		cout<<"\nMatrices Incompatible for dot porduct";
		return dots;
	}
	else
	{
		for(int i = 0; i < num_dims; i++)
		{
			if(shape[i] != a.shape[i])
			{
				cout<<"\nMatrices Incompatible for dot porduct";
				return dots;
			}
		}
		for(int i = 0; i < endval; i++)
			for(int j = 0 ; j < num_dims; j++)
			{
				if((indices[i][j] == a.indices[i][j]) && (indices[i][j+1] == a.indices[i][j+1]))
				{
					dots[i] += nz_vals[i] * a.nz_vals[i];
					break;
				}					
			}
				
	}
	return dots;
}
