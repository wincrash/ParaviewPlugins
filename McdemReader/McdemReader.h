#ifndef __McdemReader_h
#define __McdemReader_h

#include "vtkUnstructuredGridAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <map>




class McdemReader : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(McdemReader,vtkUnstructuredGridAlgorithm)
    void PrintSelf(ostream& os, vtkIndent indent);

    static McdemReader *New();

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
    McdemReader();
    ~McdemReader(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    McdemReader(const McdemReader&);  // Not implemented.
    void operator=(const McdemReader&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
