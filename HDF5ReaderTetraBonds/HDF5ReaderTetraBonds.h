#ifndef __HDF5ReaderTetraBonds_h
#define __HDF5ReaderTetraBonds_h

#include "vtkPolyDataAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <map>




class HDF5ReaderTetraBonds : public vtkPolyDataAlgorithm
{
public:
    vtkTypeMacro(HDF5ReaderTetraBonds,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    static HDF5ReaderTetraBonds *New();

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
    HDF5ReaderTetraBonds();
    ~HDF5ReaderTetraBonds(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    HDF5ReaderTetraBonds(const HDF5ReaderTetraBonds&);  // Not implemented.
    void operator=(const HDF5ReaderTetraBonds&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
