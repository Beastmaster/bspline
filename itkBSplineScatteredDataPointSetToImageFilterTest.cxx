#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <string>
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkBSplineScatteredDataPointSetToImageFilter2.h"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>



#include "bspline_function.h"
#include "SymmetricMeanOfClosestDistance.h"
//
//In this test, we approximate a 2-D scalar field.
//The scattered data is derived from a segmented 
//image.  We write the output to an image for
//comparison.
//
int itkBSplineScatteredDataPointSetToImageFilterTest1( int argc, char **argv )
{
  const unsigned int ParametricDimension = 2;
  const unsigned int DataDimension = 1;

  typedef int PixelType;
  typedef itk::Image<PixelType, ParametricDimension> InputImageType;
  typedef float RealType;
  typedef itk::Vector<RealType, DataDimension> VectorType;
  typedef itk::Image<VectorType, ParametricDimension> VectorImageType;
  typedef itk::PointSet
    <VectorImageType::PixelType, ParametricDimension> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();  

  typedef itk::ImageFileReader<InputImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  itk::ImageRegionIteratorWithIndex<InputImageType> 
    It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  
  // Iterate through the input image which consists of multivalued 
  // foreground pixels (=nonzero) and background values (=zero).
  // The foreground pixels comprise the input point set.
  
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() != itk::NumericTraits<PixelType>::Zero )
      {
      // We extract both the 2-D location of the point 
      // and the pixel value of that point.  

      PointSetType::PointType point;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );

      unsigned long i = pointSet->GetNumberOfPoints();
      pointSet->SetPoint( i, point );        

      PointSetType::PixelType V( DataDimension );
      V[0] = static_cast<RealType>( It.Get() );
      pointSet->SetPointData( i, V );
      }
    }

  
  // Instantiate the B-spline filter and set the desired parameters.
  typedef itk::BSplineScatteredDataPointSetToImageFilter2
    <PointSetType, VectorImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetSplineOrder( 3 );  
  FilterType::ArrayType ncps;  
  ncps.Fill( 4 );  
  filter->SetNumberOfControlPoints( ncps );
  filter->SetNumberOfLevels( 3 );

  // Define the parametric domain.
  filter->SetOrigin( reader->GetOutput()->GetOrigin() );
  filter->SetSpacing( reader->GetOutput()->GetSpacing() );
  filter->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );

  filter->SetInput( pointSet );

  try 
    {
    filter->Update();
    }
  catch (...) 
    {
    std::cerr << "Test 1: itkBSplineScatteredDataImageFilter exception thrown" 
              << std::endl;
    return EXIT_FAILURE;
    }
  
  // Write the output to an image.
  typedef itk::Image<RealType, ParametricDimension> RealImageType;
  RealImageType::Pointer image = RealImageType::New();
  image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  image->Allocate();
  itk::ImageRegionIteratorWithIndex<RealImageType> 
    Itt( image, image->GetLargestPossibleRegion() );
  
  for ( Itt.GoToBegin(); !Itt.IsAtEnd(); ++Itt )
    {
    Itt.Set( filter->GetOutput()->GetPixel( Itt.GetIndex() )[0] );
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS; 
};

//
//In this example, we sample a parametric curve (helix)
//and reconstruct using B-splines.
//
int itkBSplineScatteredDataPointSetToImageFilterTest2()
{
  const unsigned int ParametricDimension = 1;
  const unsigned int DataDimension = 3;

  typedef double RealType;
  typedef itk::Vector<RealType, DataDimension> VectorType;
  typedef itk::Image<VectorType, ParametricDimension> ImageType;  

  typedef itk::PointSet<VectorType, ParametricDimension> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();  

  auto sample_helix = [](std::string file_name)
  {
	  PointSetType::Pointer pointSet = PointSetType::New();
	  std::ofstream of_stream;
	  of_stream.open(file_name);

	  // Sample the helix.
	  for (RealType t = 0.0; t <= 1.0 + 1e-10; t += 0.05)
	  {
		  unsigned long i = pointSet->GetNumberOfPoints();

		  PointSetType::PointType point;
		  point[0] = t;
		  pointSet->SetPoint(i, point);

		  VectorType V;
		  V[0] = 0.25*cos(t*6.0*3.141);
		  V[1] = 0.25*sin(t*6.0*3.141);
		  V[2] = 4.0*t;

		  of_stream << V[0] << ',' << V[1] << ',' << V[2] << std::endl;

		  pointSet->SetPointData(i, V);
	  }
	  of_stream.close();
  };

  // read from file
  auto read_from_file = [&pointSet](std::string file_name)
  {

	  std::ifstream in_stream;
	  in_stream.open(file_name);
	  std::string line;
	  while (!in_stream.eof())
	  {
		  in_stream >> line;

		  auto pos = line.find(',');
		  float x = std::stof(line.substr(0,pos));
		  auto pos2 = line.find(',',pos+1);
		  float y = std::stof(line.substr(pos+1,pos2-pos));
		  float z = std::stof(line.substr(pos2+1, line.length()));

		  std::cout << "x: " << x;
		  std::cout << "\ty: " << y;
		  std::cout << "\tz: " << z << std::endl;

		  unsigned long i = pointSet->GetNumberOfPoints();
		  PointSetType::PointType point;
		  point[0] = i/100;
		  pointSet->SetPoint(i, point);

		  VectorType V;
		  V[0] = x;
		  V[1] = y;
		  V[2] = z ;

		  pointSet->SetPointData(i, V);

	  }
	  in_stream.close();
  };

  read_from_file("sample.txt");

  // Instantiate the filter and set the parameters
  typedef itk::BSplineScatteredDataPointSetToImageFilter2
     <PointSetType, ImageType>  FilterType;
  FilterType::Pointer filter = FilterType::New();
  
  // Define the parametric domain
  ImageType::SpacingType spacing;  
  spacing.Fill( 0.001 );
  ImageType::SizeType size;  
  size.Fill( static_cast<unsigned int>( 1.0/spacing[0] )+1  );
  ImageType::PointType origin;  
  origin.Fill( 0.0 );

  filter->SetSize( size );
  filter->SetOrigin( origin );
  filter->SetSpacing( spacing );
  filter->SetInput( pointSet );

  filter->SetSplineOrder( 3 );  
  FilterType::ArrayType ncps;
  ncps.Fill( 4 );  
  filter->SetNumberOfControlPoints( ncps );
  filter->SetNumberOfLevels( 20 );//5
  filter->SetGenerateOutputImage( false );

  try 
    {
    filter->Update();
    
	std::string output;

    for ( RealType t = 0.0; t <= 1.0+1e-10; t += 0.01 )
      {
		  PointSetType::PointType point;
		  point[0] = t;

		  //output.append("t:");
		  //output.append(std::to_string(t));
		  //output.append("\n");

		  VectorType V; 
		  filter->Evaluate( point, V );
		  //output.append("V:");
		  output.append(std::to_string(V[0])); output.append(" ");
		  output.append(std::to_string(V[1])); output.append(" ");
		  output.append(std::to_string(V[2]));
		  output.append("\n");

		  FilterType::GradientType G;
		  filter->EvaluateGradient( point, G );
      }
	  std::ofstream of_stream;
	  of_stream.open("spline.txt");
	  of_stream << output;
	  of_stream.close();
    }
  catch (...) 
    {
    std::cerr << "Test 2: itkBSplineScatteredDataImageFilter exception thrown" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
};


bool Test_Bspline_function()
{
	const unsigned parametericDimension = 1;
	const unsigned dataDimension = 3;

	typedef double RealType;
	typedef itk::Vector<RealType, dataDimension> VectorType;
	typedef itk::Image<VectorType, parametericDimension> ImageType;
	typedef itk::PointSet<VectorType, parametericDimension> PointSetType;

	PointSetType::Pointer pointSet = PointSetType::New();

	auto sample_helix = [&pointSet](std::string file_name)
	{
		std::ofstream of_stream;
		of_stream.open(file_name);

		// Sample the helix.
		for (RealType t = 0.0; t <= 1.0 + 1e-10; t += 0.05)
		{
			unsigned long i = pointSet->GetNumberOfPoints();

			PointSetType::PointType point;
			point[0] = t;
			pointSet->SetPoint(i, point);

			VectorType V;
			V[0] = 0.25*cos(t*6.0*3.141);
			V[1] = 0.25*sin(t*6.0*3.141);
			V[2] = 4.0*t;

			of_stream << V[0] << ',' << V[1] << ',' << V[2] << std::endl;

			pointSet->SetPointData(i, V);
		}
		of_stream.close();
	};


	return true;
}

void test3()
{
	auto bspline = new b_spline_function();

	double point[][3] = { { 0.146982, 0.202228, 0.2 },
						  { -0.0771697, 0.237792, 0.4 },
						  { -0.237723, 0.0773811, 0.6 },
						  { -0.202359, -0.146802, 0.8 },
						  { -0.000222245, -0.25, 1 },
						  { 0.202097, -0.147162, 1.2 },
						  { 0.23786, 0.0769583, 1.4 },
						  { 0.0775924, 0.237654, 1.6 },
						  { -0.146622, 0.202489, 1.8 },
						  { -0.25, 0.00044449, 2 },
						  { -0.147342, -0.201966, 2.2 },
						  { 0.0767468, -0.237928, 2.4 },
						  { 0.237585, -0.0778036, 2.6 },
						  { 0.202619, 0.146442, 2.8 },
						  { 0.000666734, 0.249999, 3 },
						  { -0.201835, 0.147521, 3.2 },
						  { -0.237997, -0.0765352, 3.4 },
						  { -0.0780148, -0.237516, 3.6 },
						  { 0.146262, -0.202749, 3.8 },
						  { 0.249998, -0.000888979, 4 } };
	std::vector<double*> ori_points;
	int num_points = 5;
	for (int i = 0; i < num_points; i++)
	{
		ori_points.push_back(point[i]);
		bspline->AddControlPoint(point[i]);
	}
	
	int p = 3;

	bspline->SetDegree(p);

	double weight[] = {1,1,4,1,1 };
	std::vector<double> weight_v(weight, weight + sizeof(weight) / sizeof(double));
	bspline->SetWeight(weight_v);

	std::vector<double> knots_v;
	auto num = b_spline_function::GenerateUniformKnots(knots_v,p,num_points);
	
	bspline->SetKnotVector(knots_v);

	std::string output;
	std::vector<double*> out_points;
	for (double t = 0.0; t <= 1.0; t += 0.005)
	{
		double* pnt = new double[3];
		//bspline->Evaluate(t,pnt);
		bspline->Evaluate_W(t,pnt);

		out_points.push_back(pnt);
		//output.append("V:");
		output.append(std::to_string(pnt[0])); output.append(" ");
		output.append(std::to_string(pnt[1])); output.append(" ");
		output.append(std::to_string(pnt[2]));
		output.append("\n");
	}
	std::ofstream of_stream;
	of_stream.open("spline.txt");
	of_stream << output;
	of_stream.close();


	//======= visualization =========//
	auto add_point = [](std::vector<double*> input, char color)
	{
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();

		for (auto it = input.begin(); it != input.end(); ++it)
		{
			points->InsertNextPoint((*it)[0],(*it)[1],(*it)[2]);
		}

		vtkSmartPointer<vtkPolyData> pointsPolydata =
			vtkSmartPointer<vtkPolyData>::New();

		pointsPolydata->SetPoints(points);

		vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
			vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
		vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
#else
		vertexFilter->SetInputData(pointsPolydata);
#endif
		vertexFilter->Update();

		vtkSmartPointer<vtkPolyData> polydata =
			vtkSmartPointer<vtkPolyData>::New();
		polydata->ShallowCopy(vertexFilter->GetOutput());

		// Setup colors
		unsigned char red[3] = { 255, 0, 0 };
		unsigned char green[3] = { 0, 255, 0 };
		unsigned char blue[3] = { 0, 0, 255 };

		vtkSmartPointer<vtkUnsignedCharArray> colors =
			vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("Colors");

		switch (color)
		{
		case 'r':
			{
				for (auto it = input.begin(); it != input.end(); ++it)
				{
					colors->InsertNextTupleValue(red);
				}
			}
			break;
		case 'b':
			{
				for (auto it = input.begin(); it != input.end(); ++it)
				{
					colors->InsertNextTupleValue(blue);
				}
			}
			break;
		case 'y':
			{
				for (auto it = input.begin(); it != input.end(); ++it)
				{
					colors->InsertNextTupleValue(green);
				}
			}
			break;
		default:
			{
				for (auto it = input.begin(); it != input.end(); ++it)
				{
					colors->InsertNextTupleValue(red);
				}
			}
			break;
		}

		polydata->GetPointData()->SetScalars(colors);
		
		return polydata;
	};


	auto dot1 = add_point(ori_points, 'r');
	auto dot2 = add_point(out_points, 'b');


	// Visualization
	// point set 1
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputData(dot1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(10);

	// point set 2
	vtkSmartPointer<vtkPolyDataMapper> mapper2 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper2->SetInputData(dot2);
#endif

	vtkSmartPointer<vtkActor> actor2 =
		vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetPointSize(5);


	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->AddActor(actor2);

	renderWindow->Render();
	renderWindowInteractor->Start();
}

void test4()
{
	double knots[] = { 0, 0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1, 1 };
	std::vector<double> knots_v(knots, knots + sizeof(knots) / sizeof(double));

	int p = 5;
	double u = 0.5;

	auto i = b_spline_function::FindSpan(p,u,knots_v);
	
	std::vector<double> xx;

	b_spline_function::BasisFunctions(i,u,p,knots_v,xx);

	return;
}


void test5()
{

	////Setup point coordinates
	//double x[3] = { 1.0, 0.0, 0.0 };
	//double y[3] = { 0.0, 1.0, 0.0 };
	//double z[3] = { 0.0, 0.0, 1.0 };
	//
	//vtkSmartPointer<vtkPoints> points =
	//	vtkSmartPointer<vtkPoints>::New();
	//points->InsertNextPoint(x);
	//points->InsertNextPoint(y);
	//points->InsertNextPoint(z);
	//
	////Create the tree
	//vtkSmartPointer<vtkKdTree> kDTree =
	//	vtkSmartPointer<vtkKdTree>::New();
	//kDTree->BuildLocatorFromPoints(points);
	//
	//double testPoint[3] = { 100.0, 0.0, 0.0 };
	//
	////Find the closest points to TestPoint
	//double closestPointDist;
	//vtkIdType id = kDTree->FindClosestPoint(testPoint, closestPointDist); //vtkKdTree::FindClosestPoint: must build locator first
	//std::cout << "The closest point is point " << id <<"     distance:  "<<closestPointDist <<std::endl;
	//

	//==== construct  polydata  0==//
	auto poly0 = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points0 =
		vtkSmartPointer<vtkPoints>::New();
	points0->InsertNextPoint(0.0, 0.0, 0.0);
	points0->InsertNextPoint(1.0, 0.0, 0.0);
	points0->InsertNextPoint(2.0, 0.0, 0.0);
	poly0->SetPoints(points0);
	
	
	//==== construct  polydata  0==//
	auto poly1 = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points1 =
		vtkSmartPointer<vtkPoints>::New();
	points1->InsertNextPoint(0.0, 0.0, 0.0);
	points1->InsertNextPoint(1.0, 0.0, 0.0);
	points1->InsertNextPoint(2.0, 0.0, 0.0);
	//points1->InsertNextPoint(0.0, 1.0, 0.0);
	//points1->InsertNextPoint(30.0, 1.0, 0.0);
	//points1->InsertNextPoint(3.0, 1.0, 0.0);
	poly1->SetPoints(points1);

	double xx = SymmetricMeanOfCloestDistance2(poly0,poly1);
	std::cout << "Distance is :  " << xx << std::endl;

}


void test_surface()
{
	double row[4][4][3] = { 
		{
			{ 1.0, 1.0, 0 },
			{ 1.0, 2.0, 1 },
			{ 1.0, 3.0, 0 },
			{ 1.0, 4.0, 0 }
		},

		{
			{ 2.0, 1.0, 0 },
			{ 2.0, 2.0, 1 },
			{ 2.0, 3.0, 0 },
			{ 2.0, 4.0, 0 }
		},

		{
			{ 3.0, 1.0, 0 },
			{ 3.0, 2.0, 1 },
			{ 3.0, 3.0, 0 },
			{ 3.0, 4.0, 0 },
		} ,
		{
			{ 4.0, 1.0, 0 },
			{ 4.0, 2.0, 1 },
			{ 4.0, 3.0, 0 },
			{ 4.0, 4.0, 0 },
		}
	};

	std::vector<std::vector<double*> > points;
	std::vector<double*>  ori_points; // for view
	for (size_t i = 0; i < 4; i++)
	{
		std::vector<double*> temp;
		for (size_t j = 0; j < 4; j++)
		{
			temp.push_back(row[i][j]);
			ori_points.push_back(row[i][j]);
		}
		points.push_back(temp);
	}

	int degrere = 2;
	std::vector<double> knots_U;
	std::vector<double> knots_V;

	b_spline_function::GenerateUniformKnots(knots_U, degrere, 4);
	b_spline_function::GenerateUniformKnots(knots_V, degrere, 4);

	auto spline = new b_spline_function();
	spline->SetDegree_U(degrere);
	spline->SetDegree_V(degrere);

	spline->SetKnotVector_U(knots_U);
	spline->SetKnotVector_V(knots_V);

	spline->AddControlPoint(points);

	std::vector<double* > surface_points;
	for (double u = 0.0; u < 1.0; u+=0.01)
	{
		for (double v = 0.0; v < 1.0; v += 0.01)
		{
			double * pnt = new double(3);
			spline->Evaluate_Surface(u,v,pnt);
			surface_points.push_back(pnt);

			//std::cout << u << "\t" << v << std::endl;
		}
	}
	std::string output;
	for (auto it = surface_points.begin(); it != surface_points.end(); ++it)
	{
		for (size_t j = 0; j < 3; j++)
		{
			output.append(std::to_string((*it)[j])); output.append(" ");
		}
		output.append("\n");
	}
	std::ofstream of_stream;
	of_stream.open("spline.txt");
	of_stream << output;
	of_stream.close();
	


	//======= visualization =========//
	auto add_point = [](std::vector<double*> input, char color)
	{
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();

		for (auto it = input.begin(); it != input.end(); ++it)
		{
			points->InsertNextPoint((*it)[0], (*it)[1], (*it)[2]);
		}

		vtkSmartPointer<vtkPolyData> pointsPolydata =
			vtkSmartPointer<vtkPolyData>::New();

		pointsPolydata->SetPoints(points);

		vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =
			vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
		vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
#else
		vertexFilter->SetInputData(pointsPolydata);
#endif
		vertexFilter->Update();

		vtkSmartPointer<vtkPolyData> polydata =
			vtkSmartPointer<vtkPolyData>::New();
		polydata->ShallowCopy(vertexFilter->GetOutput());

		// Setup colors
		unsigned char red[3] = { 255, 0, 0 };
		unsigned char green[3] = { 0, 255, 0 };
		unsigned char blue[3] = { 0, 0, 255 };

		vtkSmartPointer<vtkUnsignedCharArray> colors =
			vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("Colors");

		switch (color)
		{
		case 'r':
		{
			for (auto it = input.begin(); it != input.end(); ++it)
			{
				colors->InsertNextTupleValue(red);
			}
		}
		break;
		case 'b':
		{
			for (auto it = input.begin(); it != input.end(); ++it)
			{
				colors->InsertNextTupleValue(blue);
			}
		}
		break;
		case 'y':
		{
			for (auto it = input.begin(); it != input.end(); ++it)
			{
				colors->InsertNextTupleValue(green);
			}
		}
		break;
		default:
		{
			for (auto it = input.begin(); it != input.end(); ++it)
			{
				colors->InsertNextTupleValue(red);
			}
		}
		break;
		}

		polydata->GetPointData()->SetScalars(colors);

		return polydata;
	};


	auto dot1 = add_point(ori_points, 'r');
	auto dot2 = add_point(surface_points, 'b');


	// Visualization
	// point set 1
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper->SetInputData(dot1);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(10);

	// point set 2
	vtkSmartPointer<vtkPolyDataMapper> mapper2 =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
#else
	mapper2->SetInputData(dot2);
#endif

	vtkSmartPointer<vtkActor> actor2 =
		vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
	actor2->GetProperty()->SetPointSize(5);


	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->AddActor(actor2);

	renderWindow->Render();
	renderWindowInteractor->Start();

}

int main( int argc, char **argv )
{
	test_surface();
	//test5();
	//test4();
	//test3();
    //bool test1 = itkBSplineScatteredDataPointSetToImageFilterTest1( argc, argv );
    //bool test2 = itkBSplineScatteredDataPointSetToImageFilterTest2();
  
	return 0;
}







