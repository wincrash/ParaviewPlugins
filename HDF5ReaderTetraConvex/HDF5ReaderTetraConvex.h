#ifndef __HDF5ReaderTetraConvex_h
#define __HDF5ReaderTetraConvex_h

#include "vtkUnstructuredGridAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <map>




class HDF5ReaderTetraConvex : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(HDF5ReaderTetraConvex,vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    static HDF5ReaderTetraConvex *New();

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
    HDF5ReaderTetraConvex();
    ~HDF5ReaderTetraConvex(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    HDF5ReaderTetraConvex(const HDF5ReaderTetraConvex&);  // Not implemented.
    void operator=(const HDF5ReaderTetraConvex&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
