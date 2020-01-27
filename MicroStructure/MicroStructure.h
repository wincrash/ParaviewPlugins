#ifndef __MicroStructure_h
#define __MicroStructure_h

#include "vtkSmartPointer.h" // compiler errors if this is forward declared


#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkExecutive.h"

class vtkUnstructuredGrid;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;


class MicroStructure : public vtkUnstructuredGridAlgorithm
{
 public:
    double CELL_SIZE;
    vtkSetMacro(CELL_SIZE, double);
    double RADIUS;
    vtkSetMacro(RADIUS, double);

  static MicroStructure *New();
  vtkTypeMacro(MicroStructure, vtkUnstructuredGridAlgorithm);
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
  MicroStructure();
  ~MicroStructure();

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );

  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

 private:
};
#endif
