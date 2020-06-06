#ifndef __HDF5ReaderPeridynamics_h
#define __HDF5ReaderPeridynamics_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyDataAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <map>




class HDF5ReaderPeridynamics : public vtkPolyDataAlgorithm
{
public:
    vtkTypeMacro(HDF5ReaderPeridynamics,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    static HDF5ReaderPeridynamics *New();

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
    HDF5ReaderPeridynamics();
    ~HDF5ReaderPeridynamics(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    HDF5ReaderPeridynamics(const HDF5ReaderPeridynamics&);  // Not implemented.
    void operator=(const HDF5ReaderPeridynamics&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
