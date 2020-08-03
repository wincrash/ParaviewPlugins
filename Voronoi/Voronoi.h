#ifndef __Voronoi_h
#define __Voronoi_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared

#include "vtkUnstructuredGridAlgorithm.h" //superclass
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkExecutive.h"

class vtkUnstructuredGrid;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;



class Voronoi : public vtkUnstructuredGridAlgorithm
{
 public:
  static Voronoi *New();
  vtkTypeMacro(Voronoi, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Specify the source object. This is the object that will be moved during the transformation.
//  vtkUnstructuredGrid *GetSource();

  // Description:
  // Specify the target object. This is the object that will stay in place.
//  vtkUnstructuredGrid *GetTarget();

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
