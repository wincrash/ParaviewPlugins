#ifndef __McsdemReader_h
#define __McsdemReader_h

#include "vtkUnstructuredGridAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <map>




class McsdemReader : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(McsdemReader,vtkUnstructuredGridAlgorithm)
    void PrintSelf(ostream& os, vtkIndent indent);

    static McsdemReader *New();

    // Description:
    // Specify file name of the .abc file.
    vtkSetStringMacro(FileName)
    vtkGetStringMacro(FileName)

    void  SetDirectoryName(const char* dn)
    {
        this->DirectoryName = strdup(vtksys::SystemTools::GetParentDirectory(dn).c_str());
    }
    vtkGetStringMacro(DirectoryName)
protected:
    McsdemReader();
    ~McsdemReader(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    McsdemReader(const McsdemReader&);  // Not implemented.
    void operator=(const McsdemReader&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
