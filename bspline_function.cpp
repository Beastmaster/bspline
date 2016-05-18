/*
Author: QIN Shuo
Date:   2016/5/16

Basic spline function

*/



#include "bspline_function.h"




b_spline_function::b_spline_function()
{
	m_P_degree = 3;
}

b_spline_function::~b_spline_function()
{
	for (auto it = m_Control_Points.begin(); it != m_Control_Points.end(); ++it)
	{
		delete *it;
	}
}




void b_spline_function::AddControlPoint(double* point)
{
	double* pp = new double[3];
	memcpy(pp,point,3*sizeof(double));

	this->m_Control_Points.push_back(pp);
}

void b_spline_function::SetKnotVector(std::vector<double> in)
{
	m_knot_vector = in;
}

void b_spline_function::SetDegree(int p)
{
	this->m_P_degree = p;
}


void b_spline_function::Evaluate(double u ,double * coor)
{
	std::vector<double> coeff;
	auto span = FindSpan(m_P_degree,u,m_knot_vector);
	BasisFunctions(span,u,m_P_degree,m_knot_vector,coeff);

	memset(coor,0,3*sizeof(double));

	for (size_t i = 0; i <= m_P_degree; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			// n-1 basic functions
			coor[j] = coor[j] + coeff[i] * m_Control_Points[span - m_P_degree + i ][j];
			//std::cout << i << "   " << j <<"   "<< span << std::endl;
		}
	}

}



void b_spline_function::AddControlPoint(std::vector<std::vector<double*> > input)
{
	m_Sur_ControlPoints = input;
}
void b_spline_function::SetDegree_U(int p)
{
	m_P_U = p;
}
void b_spline_function::SetDegree_V(int q)
{
	m_P_V = q;
}

void b_spline_function::SetKnotVector_U(std::vector<double> in)
{
	m_knot_U = in;
}
void b_spline_function::SetKnotVector_V(std::vector<double> in)
{
	m_knot_V = in;
}

void b_spline_function::Evaluate_Surface(double u, double v, double* out)
{
	std::vector<double> coeff_U;
	std::vector<double> coeff_V;
	int uspan = this->FindSpan(m_P_U, u, m_knot_U);
	int vspan = this->FindSpan(m_P_V, v, m_knot_V);

	BasisFunctions(uspan, u, m_P_U, m_knot_U, coeff_U);
	BasisFunctions(vspan, v, m_P_V, m_knot_V, coeff_V);
	
	int uid = uspan - m_P_U;
	
	//auto S = new double(3);
	memset(out, 0, 3 * sizeof(double));

	for (size_t l = 0; l <= m_P_V;l++)
	{
		double* temp = new double(3); memset(temp,0,3* sizeof(double));
		int vid = vspan - m_P_V + l;
		for (size_t k = 0; k <= m_P_U; k++)
		{
			for (size_t i = 0; i < 3; i++)
				temp[i] = temp[i] + coeff_U[k] * m_Sur_ControlPoints[uid + k][vid][i];
		}

		for (size_t i = 0; i < 3; i++)
		{
			out[i] = out[i] + coeff_V[l] * temp[i];
		}
	}

}







/*
Weight curve
*/
void b_spline_function::SetWeight(std::vector<double> weight)
{
	m_Weight = weight;
}

/*
Allocate memory for out before call
*/
void b_spline_function::Evaluate_W(double u, double* out)
{
	std::vector<double> coeff;
	auto span = FindSpan(m_P_degree, u, m_knot_vector);
	BasisFunctions(span, u, m_P_degree, m_knot_vector, coeff);

	double* Cw = new double(3);
	memset(Cw, 0, 3 * sizeof(double));
	memset(out, 0, 3 * sizeof(double));

	double w = 0.0;
	for (size_t i = 0; i <= m_P_degree; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			// n-1 basic functions
			Cw[j] = Cw[j] + coeff[i] * m_Control_Points[span - m_P_degree + i][j] * m_Weight[span - m_P_degree + i];
			//std::cout << i << "   " << j <<"   "<< span << std::endl;
		}
		w = w + coeff[i] * m_Weight[span - m_P_degree + i];
	}
	for (int j = 0; j < 3; j++)
	{
		out[j] = Cw[j]/w;
	}
}


/*
m = n + p + 1
num_Knots = m+1

n: number of control points -1
p: degree

u range (0,1)

return knots vector length
*/
int b_spline_function::GenerateUniformKnots(std::vector<double> &knots,int degree, int num_controlPoints)
{
	int n = num_controlPoints -1 ; //count from 0
	int p = degree;
	int m = n + p + 1;

	knots.clear();
	for (size_t i = 0; i < p+1; i++)
	{
		knots.push_back(0.0);
	}

	double u0 = 1.0 / (m-2*p);
	for (size_t i = p+1; i < m-p ; i++)
	{
		knots.push_back((i - p)*u0);
	}

	for (size_t i = m-p; i < m+1; i++)
	{
		knots.push_back(1.0);
	}


	//test
	for (auto i = knots.begin(); i != knots.end(); ++i)
	{
		std::cout<< *i<<" , ";
	}std::cout << std::endl;

	return m + 1; //return number of knots 
}



