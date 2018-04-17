/*

 #     # ###  #####  ######     #    ######  ####### ######  ####### #     #
 #     #  #  #     # #     #   # #   #     #    #    #     # #       ##   ##
 #     #  #  #       #     #  #   #  #     #    #    #     # #       # # # #
 #     #  #   #####  ######  #     # ######     #    #     # #####   #  #  #
 #   #   #        # #       ####### #   #      #    #     # #       #     #
 # #    #  #     # #       #     # #    #     #    #     # #       #     #
 #    ###  #####  #       #     # #     #    #    ######  ####### #     #


 * Author       : Ruslan Pacevic
 * E-mail       : rpa@sc.vgtu.lt
 * Vendor       : VGTU
 * Home page    : http://lsl.vgtu.lt/vispartdem
 */
#ifndef vtkCellBased3Dver2_H_
#define vtkCellBased3Dver2_H_

#include "vtkPolyDataAlgorithm.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include <sstream>
#include "vtkCellData.h"
#include <map>
#include <algorithm>
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPyramid.h"
#include "vtkMath.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkExtractCells.h"
#include "Timeris.h"
using namespace std;
namespace vispartdem {
namespace cellbased3DRev2 {
typedef int INT;
typedef pair<INT, INT> _line;
typedef pair<_line, INT> _trik;
typedef pair<_trik, INT> _pyramid;
typedef pair<pair<pair<pair<pair<INT, INT> , INT> , INT> , INT> , INT> _keturkampis;

typedef struct MyStructure {
        vector<INT> ids;
        vector<INT> cellIDs;

        MyStructure() {
                ids.reserve(6);
                cellIDs.reserve(12);
        }
        MyStructure(const MyStructure &c) {
                ids.resize(c.ids.size(), -1);
                cellIDs.resize(c.cellIDs.size(), -1);
                copy(c.ids.begin(), c.ids.end(), ids.begin());
                copy(c.cellIDs.begin(), c.cellIDs.end(), cellIDs.begin());
        }
} MyStruct;

class VTK_GRAPHICS_EXPORT vtkCellBased3Dver2: public vtkUnstructuredGridAlgorithm {
public:
        static vtkCellBased3Dver2 *New();vtkTypeRevisionMacro(vtkCellBased3Dver2, vtkUnstructuredGridAlgorithm)
        ;
        void PrintSelf(ostream& os, vtkIndent indent);
        void SetStateArray(std::string ResultArrayName);
        const char *getStateArray();
        //vtkSetStringMacro(StateArray);
        //vtkGetStringMacro(StateArray);
        void printArray(string name,double* arr,int l)
        {
                cout<<name.c_str()<<" ";
                for(int k=0;k<l;k++)
                {
                        cout<<arr[k]<<" ";
                }
                cout<<endl;
        }


        void GetLineID(_line &line, INT id1, INT id2) {
                if (id1 < id2) {
                        line.first = id1;
                        line.second = id2;

                } else {
                        line.first = id2;
                        line.second = id1;
                }
        }

        void GetTriangleID(_trik &trike, INT id1, INT id2, INT id3) {
                int a[3];
                a[0] = id1;
                a[1] = id2;
                a[2] = id3;
                sort(a, a + 3);
                trike.first.first = a[0];
                trike.first.second = a[1];
                trike.second = a[2];
        }
        void GetPyramidID(_pyramid &piramid, INT id1, INT id2, INT id3, INT id4) {
                int a[4];
                a[0] = id1;
                a[1] = id2;
                a[2] = id3;
                a[3] = id4;
                sort(a, a + 4);
                piramid.first.first.first = a[0];
                piramid.first.first.second = a[1];
                piramid.first.second = a[2];
                piramid.second = a[3];
        }

        void GetKeturkampis(_keturkampis &ket, INT id1, INT id2, INT id3, INT id4, INT id5, INT id6) {
                int a[6];
                a[0] = id1;
                a[1] = id2;
                a[2] = id3;
                a[3] = id4;
                a[4] = id5;
                a[5] = id6;
                sort(a, a + 6);
                ket.first.first.first.first.first = a[0];
                ket.first.first.first.first.second = a[1];
                ket.first.first.first.second = a[2];
                ket.first.first.second = a[3];
                ket.first.second = a[4];
                ket.second = a[5];
        }

protected:
        vtkCellBased3Dver2();

        ~vtkCellBased3Dver2() {
        }

        // Usual data generation method
        virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
        virtual int FillInputPortInformation(int port, vtkInformation *info);

private:
        vtkCellBased3Dver2(const vtkCellBased3Dver2&); // Not implemented.
        void operator=(const vtkCellBased3Dver2&); // Not implemented.
        string StateArray;
        vtkCellArray*piramid;
        bool onetime;
        vector<MyStruct> strukturos;
    	vtkUnstructuredGrid*cells;

};
}
}
#endif
