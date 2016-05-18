/*
Author: QIN Shuo
Date:   2016/5/16

Basic spline function


[0,p] [p+1,m-p-1] [m-p,m]
p+1    m-2p-1       p+1

*/

#include <vector>
#include <iostream>

class b_spline_function
{
public:
	b_spline_function();
	~b_spline_function();

#pragma region Essential_Functions_Curve

	void AddControlPoint(double*);
	/* 3 degree by default */
	void SetDegree(int p=3);
	void SetKnotVector(std::vector<double> in);
	/*
	Evaluate coordinate in u
	*/
	void Evaluate(double u , double* out);
#pragma endregion

#pragma region Essential_Functions_Curve
	void AddControlPoint(std::vector<std::vector<double*> >);

	void SetDegree_U(int p);
	void SetDegree_V(int q);

	void SetKnotVector_U(std::vector<double> in);
	void SetKnotVector_V(std::vector<double> in);
	
	void Evaluate_Surface(double u, double v, double* out);
#pragma endregion


#pragma region Optional_Functions 
	// Set weight vector, lengh must match points
	// Have no test yet
	void SetWeight(std::vector<double> weight);
	void Evaluate_W(double u, double* out);
#pragma endregion


#pragma region GeneralFunctions
	// generate a uniform knots list
	static int GenerateUniformKnots(std::vector<double> &knots,int degree, int num_controlPoints);
#pragma endregion


	/*
	Knot vector
	U:  {a,...,a,u(p+1),...,u(m-p-1),b,...,b}
		 -------						 -------
			p+1							p+1
	p: p-th degree b-spline function
	m+1: number of knots, including first (p+1) a and last (p+1) b
	n+1 basis functions, span index
	n = m-p-1 = length - p - 2 , number of control points

	*/
	static int FindSpan( int p, double u, std::vector< double > U)
	{
		int n = U.size() - p - 2;
		if (u == U[n+1])
		{
			return n;
		}
		int low = p;
		int high = n+1 ;  // ATTENTION!!: different from the algorithm from the NURBS book
		int mid = (low + high) / 2;
		//while (u<U[mid] || u>U[mid+1])
		while (mid != low)
		{
			if (u < U[mid])
				high = mid;
			else
				low = mid;
			mid = (low + high) / 2;
		}
		return mid;
	};

	/*
	i: knot span
	u: value
	p: degree of spline
	U: knot vector
	N: 
	*/
	static void BasisFunctions(int i,double u, int p, std::vector<double> U, std::vector<double> &N)
	{
		//auto i = FindSpan(p, u, U);

		N.clear();
		N.push_back(1.0);

		// left[j]  = u - u(i+1-j)
		// right[j] = u(i+j) - u
		std::vector< double > left;
		std::vector< double > right;

		for (size_t j = 0; j <= p; j++)
		{
			left.push_back(u - U[i + 1 - j]);
			right.push_back(U[i + j] - u);
		}

		for (int j = 1; j <= p;j++)
		{
			double saved = 0.0;
			for (size_t r = 0; r < j; r++)
			{
				auto temp = N[r] / (right[r + 1] + left[j - r]);
				N[r] = saved + right[r + 1] * temp;
				saved = left[j - r] * temp;
			}
			N.push_back(saved);
		}
	};



private:
	//====  bspline curve ====//
	std::vector< double > m_knot_vector;
	std::vector< double*> m_Control_Points;
	int					  m_P_degree;
	std::vector< double > m_Weight;

	//===== bspline surface =====//
	std::vector < std::vector<double*> > m_Sur_ControlPoints;
	std::vector< double > m_knot_U;
	std::vector< double > m_knot_V;
	int					  m_P_U;//degree
	int					  m_P_V;//degree


protected:


};














