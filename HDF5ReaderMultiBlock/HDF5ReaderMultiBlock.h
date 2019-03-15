#ifndef __HDF5ReaderMultiBlock_h
#define __HDF5ReaderMultiBlock_h

#include "vtkPolyDataAlgorithm.h"
#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>
#include <map>
#include "vtkMultiBlockDataSetAlgorithm.h"


class HDF5ReaderMultiBlock : public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeMacro(HDF5ReaderMultiBlock,vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    static HDF5ReaderMultiBlock *New();

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
    HDF5ReaderMultiBlock();
    ~HDF5ReaderMultiBlock(){}

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);


private:
    vtkIdType findClosestSolutionIndex(double Time);
    HDF5ReaderMultiBlock(const HDF5ReaderMultiBlock&);  // Not implemented.
    void operator=(const HDF5ReaderMultiBlock&);  // Not implemented.

    char* FileName;
    char* DirectoryName;

    vtkSmartPointer<vtkDoubleArray> times;
    vtkSmartPointer<vtkVariantArray> timesNames;
    std::vector<std::string> filenames;


};

#endif
