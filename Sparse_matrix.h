#include<iostream>
#include<math.h>

using namespace std;

class Sparse_matrix
{

private:
//The data types are specified int initially for simplicity. These will be converted to Template format soon. 
//Values in the sparse matrix are considered to be real.

	int num_dims;               // No. of dimensions or Dimensionality
	int *shape;                 // Shape of the matrix. For e.g. 3 X 2 X 2
	int num_nzvals;             // No. of non-zero values in the sparse matrix
	double *nz_vals;            // Array of non-zero values
	int **ptr;                  // Pointer to the values(compressed format)
	int **indices;  	    // Indices for non- zero values

public:
// Member functions
	void printdata();
//	void sortindices();                         // Sorting indices based on a particular axis
	void sortindices(int axis);                 // Sorting indices - row first format

	// Arithmetic Operations
	double* cumprod();	     // To be changed as a single function that may take axis as arguement
	double* cumsum();	     // To be changed as a single function that may take axis as arguement
	double max();				    // Maximum value of the matrix
	int** argmax(int axis);                // Index of the maximum value, in the whole matrix or an axis
	double min();				    // Minimum value in the matrix
	int* argmin();                              // Index of the minimum value
	double mean();                              // Mean value of the matrix
	void sort();                                // Sorted matrix, in-place sorting
//	int** argsort();
	void __neg__();                             // -matrix
	void __pos__();                             // +matrix
	void __abs__();                             // |matrix|
//      Sparse_matrix square();
        Sparse_matrix sq_root();                    // Square root of each element
        void square();                              // matrix^2
        Sparse_matrix addition(Sparse_matrix);      // Matrix Addition
        double* dotprod(Sparse_matrix &);           // Dot product
	Sparse_matrix slice(int *);
	//Operation on values
	void clip(double min, double max);          // Restrict the values of the matrix between the min and max values
	void swapaxes(int axis1, int axis2);        // Swap the values in two axes
	void flatten();                             // Convert the matrix into an array
	double getitem(int *index);                 // Get the value corresponding to the index
	void setitem(int *index, double value);     // Set the value of an item corresponding to the index
	void getshape();                            // Get the shape of the matrix
	void addval(double value,int *index);       // add value at the particular index(same as setitem)
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
	sortindices(-1);
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

/*void Sparse_matrix :: sortindices(int axis)
{
	int k = axis, tempval=0;
	int *temp;
	temp = new int[num_dims];
	for(int i = 0; i < num_nzvals; i++)
                        for(int j = i+1 ; j < num_nzvals; j++)
                        {
                                if((indices[i][k+1] == indices[j][k+1]) && (indices[i][k] > indices[j][k]))
                                {
	                                temp = indices[i];
                                        indices[i] = indices[j];
                                        indices[j] = temp;

                                        tempval = nz_vals[i];
                                        nz_vals[i] = nz_vals[j];
                                        nz_vals[j] = tempval;
                                }
                         
			}
}
*/
void Sparse_matrix :: sortindices(int axis = -1)
{
	if(axis == -1)
	{
	int k = num_dims - 1, tempval = 0;
	int *temp;
	temp = new int[num_dims];
	while(k >= 1)
	{
		for(int i = 0; i < num_nzvals; i++)
		for(int j = i+1 ; j < num_nzvals; j++)
		{
			if(k == num_dims - 1)
			{
				if(indices[i][k] > indices[j][k])
				{	
					temp = indices[i];
					indices[i] = indices[j];
					indices[j] = temp;
					
					tempval = nz_vals[i];
					nz_vals[i] = nz_vals[j];
					nz_vals[j] = tempval;
				}
			}
			else
			{
				if(k == 1)
				{
					if(indices[i][k+1] == indices[j][k+1])
					{	if(indices[i][k-1] > indices[j][k-1])
						{
							temp = indices[i];
							indices[i] = indices[j];
							indices[j] = temp;

							tempval = nz_vals[i];
       	        	                                nz_vals[i] = nz_vals[j];
        	                                        nz_vals[j] = tempval;
						}
						if(indices[i][k-1] == indices[j][k-1])
						{
							if(indices[i][k] > indices[j][k])
							{
								temp = indices[i];
								indices[i] = indices[j];
								indices[j] = temp;
								
								tempval = nz_vals[i];
		                                                nz_vals[i] = nz_vals[j];
                       			                        nz_vals[j] = tempval;
							}
						}
					}
				}
				else
				{
					if((indices[i][k+1] == indices[j][k+1]) && (indices[i][k] > indices[j][k]))
					{	
						temp = indices[i];
						indices[i] = indices[j];
						indices[j] = temp;
				
						tempval = nz_vals[i];
                                       	        nz_vals[i] = nz_vals[j];
                                               	nz_vals[j] = tempval;
					}
				}
			}
		}
		k--;	
	}
	}
	else
	{
	int k = axis, tempval=0;
        int *temp;
        temp = new int[num_dims];
        for(int i = 0; i < num_nzvals; i++)
                        for(int j = i+1 ; j < num_nzvals; j++)
                        {
                                if((indices[i][k+1] == indices[j][k+1]) && (indices[i][k] > indices[j][k]))
                                {
                                        temp = indices[i];
                                        indices[i] = indices[j];
                                        indices[j] = temp;

                                        tempval = nz_vals[i];
                                        nz_vals[i] = nz_vals[j];
                                        nz_vals[j] = tempval;
                                }

                        }
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
	
double*  Sparse_matrix :: cumprod()	// Product of all non-zero elements - Cumulative Product without axis 
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
/*
double* Sparse_matrix :: max(int axis)
{
	double *max;
	int vals = 1;
	for(int i =0; i < numdims; i++)
	{
		if(i != axis)
			vals = vals * shape[i];
	}
	max = new double[vals];
	int k = 0;
	max[0] = nz_vals[0];
	for(int i = 0; i < num_nzvals; i++)
	{
		
		
	}
	return max;
}
*/
/*
int** Sparse_matrix :: argmax(int axis = -1)          // Returns the indices of the maximum value
{
	int **argsmax;
	sortincides(axis);
	int toalloc = 1;
	if(axis != -1)
	{
		for(int i = 0; i < num_nzvals; i++)
		{
			if(indices[i][axis] != indices[i+1][axis])
				toalloc += 1;
			else
				continue;
		}
		argsmax = new int*[toalloc];
		for(int i = 0; i < toalloc; i++)
			argsmax[i] = new int[num_dims];
		int k = 0,i = 0;
		double max;
		for(int i = 0; i < num_nzvals; i++)
		{
			max = nz_vals[i];
			if(axis == 0)
			{
			}
			else
			{	
				if(indices[i][axis+1] != indices[i+1][axis+1])
				{
					argsmax[k] = indices[i+1];
					k++;
				}
				else
				{
					if(indices[i][axis] != indices[i+1][axis])
					{
						argsmax[k] = indices[i];
						k++;
					}
					else
					{
						if(nz_vals[i] > max)
							argsmax[k] = indices[i];
							k++;
					}
				}
			}
		}
		return argsmax;
	}
	else
	{
		double max = nz_vals[0];
        	for(int i = 1; i < num_nzvals; i++)
        	{
                	if(nz_vals[i] > max)
			{
                        	max = nz_vals[i];
				*argsmax = indices[i];
			}
        	}
	}
        return argsmax;
}
*/
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

Sparse_matrix Sparse_matrix :: addition(Sparse_matrix b)	//Matrix addition - index-wise
{
	int flag = 0;
	// Step 1: Compatibility Check
	if(num_dims != b.num_dims) 
	{
		cout<<"\n Matrices Incompatible for addition";
		return b;
	}
	else
	{ // This block will be removed once != operator is overloaded
		for(int i = 0; i < num_dims; i++)
		{
			if(shape[i] != b.shape[i])
			{
				cout<<"\nMatrices Incompatible for addition";
				return b;
			}
		}
		for(int i = 0; i < b.num_nzvals; i++)
		{
			for(int j = 0; j < num_nzvals; j++)
			{
				flag = 0;
				for(int k =0; k < num_dims; k++)
				{
					if(b.indices[i][k] != indices[j][k])
					{
						flag = 1; break;
					}
				}
				if(flag != 1)
				{
					b.nz_vals[i] += nz_vals[j];
					nz_vals[j] = 0;
				}
			}
		}
		for(int i =0; i < num_nzvals; i++)
		{
			if(nz_vals[i] != 0)
			{
				b.num_nzvals++;
				b.nz_vals[b.num_nzvals] = nz_vals[i];
				b.indices[b.num_nzvals] = indices[i];
			}
			
		}
	}
	return b;
}


double* Sparse_matrix :: dotprod(Sparse_matrix &a) 
{
	double *dots;
	int endval, l=0, flag =0;
	if(num_dims != a.num_dims)
	{
		cout<<"\nMatrices Incompatible for dot product";
		return dots;
	}
	else
	{
		for(int i = 0; i < num_dims; i++)
		{
			if(shape[i] != a.shape[i])
			{
				cout<<"\nMatrices Incompatible for dot product";
				return dots;
			}
		}
		if(num_nzvals < a.num_nzvals)
			endval = num_nzvals;
		else
			endval = a.num_nzvals;
		dots = new double[endval];
		for(int i = 0; i < endval; i++)
                	dots[i] = 0;
		for(int i = 0; i < num_nzvals; i++)
		{
			for(int j = 0 ; j < a.num_nzvals; j++)
			{
				flag = 0;
				for(int k = 0; k < num_dims; k++)
					if(indices[i][k] != a.indices[j][k])
					{
						flag = 1;break;
					}
				if(flag != 1)
				{
					dots[l] += nz_vals[i] * a.nz_vals[j];
					l++;
					break;
				}				
			}
		}

	}
	return dots;
}

Sparse_matrix Sparse_matrix :: slice(int *axes)
{
	Sparse_matrix a(num_dims, num_nzvals,nz_vals,shape,indices)
	for(int i = 0;i < num_nzvals; i++)
	{
		for(int j = 0; j < num_dims; j++)
		{
			
		}
	}
	return a;
}
		
