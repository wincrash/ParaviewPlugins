#ifndef __Voronoi_h
#define __Voronoi_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkPolyDataAlgorithm.h"
#include "vtkExecutive.h"

class vtkPolyData;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;

#include "vtkMath.h"

class Voronoi : public vtkPolyDataAlgorithm
{
 public:
  static Voronoi *New();
  vtkTypeMacro(Voronoi, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Specify the source object. This is the object that will be moved during the transformation.
//  vtkPolyData *GetSource();

  // Description:
  // Specify the target object. This is the object that will stay in place.
//  vtkPolyData *GetTarget();

  void AddSourceConnection(vtkAlgorithmOutput* input);
  void RemoveAllSources();

 protected:
  Voronoi();
  ~Voronoi();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

 private:
};
#endif
