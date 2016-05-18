/*
* File:		main.cxx
* Description:	A utility to build a graph from tirangled surface, and to extract connected points.
* Author:	Zhang Teng, PhD Candidate
* Organization:	Department of Imaging and Interventional Radiology, Chinese University of Hong Kong
* Mailbox:	zhangteng630@gmail.com
* License:	GNU GPL
*
* Created on January 15, 2016, 10:21 AM
*/

#ifndef _SYMMETRICMEANOFCLOESTDISTANCE_H_
#define _SYMMETRICMEANOFCLOESTDISTANCE_H_


#include <algorithm>
#include <functional>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkKdTree.h>
#include <vtkVertexGlyphFilter.h>



double SymmetricMeanOfCloestDistance(
	vtkSmartPointer<vtkPolyData> poly0,
	vtkSmartPointer<vtkPolyData> poly1
	)
{
	//double calculate = [](vtkSmartPointer<vtkPoints> poly0, vtkSmartPointer<vtkPoints> poly1)
	//{
		//Create the tree
		vtkSmartPointer<vtkKdTree> kDTree =
			vtkSmartPointer<vtkKdTree>::New();
		kDTree->BuildLocatorFromPoints(poly0->GetPoints());

		double return_distance = 0.0;
		// Write all of the coordinates of the points in the vtkPolyData to the console.
		for (vtkIdType i = 0; i < poly1->GetNumberOfPoints(); i++)
		{
			double p[3];
			poly1->GetPoint(i, p);
		
			//Find the closest points to TestPoint
			double closestPointDist;
			vtkIdType id = kDTree->FindClosestPoint(p, closestPointDist);
			return_distance += closestPointDist;

			std::cout << "p: " << p[0] << "  " << p[1] << "  " << p[2] << "  closest:  " << closestPointDist << std::endl;
		}
		double mean_distance1 = return_distance / poly1->GetNumberOfPoints();
		std::cout << "man:  " << mean_distance1 << std::endl;
	//};


		//=========================================//
		//Create the tree
		vtkSmartPointer<vtkKdTree> kDTree2 =
			vtkSmartPointer<vtkKdTree>::New();
		kDTree2->BuildLocatorFromPoints(poly1->GetPoints());
		double return_distance2 = 0.0;
		// Write all of the coordinates of the points in the vtkPolyData to the console.
		for (vtkIdType i = 0; i < poly0->GetNumberOfPoints(); i++)
		{
			double p[3];
			poly0->GetPoint(i, p);

			//Find the closest points to TestPoint
			double closestPointDist;
			vtkIdType id = kDTree2->FindClosestPoint(p, closestPointDist);

			return_distance2 += closestPointDist;

			std::cout <<"p: " << p[0] << "  " << p[1] << "  " << p[2] <<"  closest:  " << closestPointDist<<std::endl;
		}
		double mean_distance2 = return_distance2 / poly0->GetNumberOfPoints();
		std::cout << "mean:  " << mean_distance2 << std::endl;

		//return
		return (mean_distance1 + mean_distance2) / 2;
}

double SymmetricMeanOfCloestDistance2(
	vtkSmartPointer<vtkPolyData> poly0,
	vtkSmartPointer<vtkPolyData> poly1
	)
{
	// calculate distance: points in poly0 --->>  poly1
	auto calculate = [](vtkSmartPointer<vtkPolyData> poly0, vtkSmartPointer<vtkPolyData> poly1)
	{
		auto min_distance = [](double* coordinate, vtkSmartPointer<vtkPoints> points)
		{
			double distance = DBL_MAX;
			for (size_t i = 0; i < points->GetNumberOfPoints(); i++)
			{
				double p[3];
				points->GetPoint(i,p);

				auto cal_point_distance = [](double* p1, double* p2)
				{
					double dis = 0.0;
					for (size_t i = 0; i < 3; i++)
					{
						dis += (p1[i] - p2[i]) *(p1[i] - p2[i]);
					}
					dis = sqrt(dis);
					return dis;
				};

				auto temp = cal_point_distance(coordinate, p);
				if (distance >= temp)
				{
					distance = temp;
				}
			}
			return distance;
		};

		double sum = 0.0;
		for (size_t i = 0; i < poly0->GetNumberOfPoints(); i++)
		{
			double p[3];
			poly0->GetPoint(i,p);

			double min_ = min_distance(p,poly1->GetPoints());
			sum += min_;
		}

		return sum / poly0->GetNumberOfPoints();
	};


	double distance1 = calculate(poly0,poly1);
	double distance2 = calculate(poly1,poly0);

	return (distance1 + distance2) / 2;
}

#endif