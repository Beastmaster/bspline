cmake_minimum_required(VERSION 2.8)

set(PRJ itk_bspline)

project(${PRJ})
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
  
 set(
 SRCS
 itkBSplineScatteredDataPointSetToImageFilterTest.cxx
 itkBSplineScatteredDataPointSetToImageFilter2.h
 itkBSplineScatteredDataPointSetToImageFilter2.txx
 
 #itkPointSetToImageFilter2.h
 #itkPointSetToImageFilter2.txx
 #
 #itkBSplineKernelFunction2.h
 #itkBSplineKernelFunction2.txx
 #
 #itkKernelFunction.h
 
 )
 
 
add_executable(${PRJ} MACOSX_BUNDLE ${SRCS})

target_link_libraries(${PRJ} ITKReview  ${ITK_LIBRARIES})
target_link_libraries(${PRJ} ${VTK_LIBRARIES})



