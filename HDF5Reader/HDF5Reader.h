#ifndef __HDF5Reader_h
#define __HDF5Reader_h
 
#include "vtkPolyDataAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <map>

 
class HDF5Reader : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(HDF5Reader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static HDF5Reader *New();
 
  // Description:
  // Specify file name of the .abc file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

void  SetDirectoryName(const char* dn)
  {
    this->DirectoryName = strdup(vtksys::SystemTools::GetParentDirectory(dn).c_str());
  }
   vtkGetStringMacro(DirectoryName);
protected:
  HDF5Reader();
  ~HDF5Reader(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);

 
private:
  vtkIdType findClosestSolutionIndex(double Time);
  HDF5Reader(const HDF5Reader&);  // Not implemented.
  void operator=(const HDF5Reader&);  // Not implemented.
 
  char* FileName;
  char* DirectoryName;

  vtkSmartPointer<vtkDoubleArray> times;
  vtkSmartPointer<vtkVariantArray> timesNames;


};
 
#endif
