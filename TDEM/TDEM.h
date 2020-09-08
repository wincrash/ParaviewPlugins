#ifndef __TDEM_h
#define __TDEM_h

#include "vtkUnstructuredGridAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkVariantArray.h>
#include <map>




class TDEM : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(TDEM,vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    static TDEM *New();

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
    TDEM();
    ~TDEM(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    TDEM(const TDEM&);  // Not implemented.
    void operator=(const TDEM&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
