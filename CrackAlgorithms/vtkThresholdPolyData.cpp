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
#include "vtkThresholdPolyData.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"


vtkStandardNewMacro(vtkThresholdPolyData);


vtkThresholdPolyData::vtkThresholdPolyData() {
    this->LowerThreshold = 0.0;
    this->UpperThreshold = 1.0;
    this->AllScalars = 1;
    this->AttributeMode = -1;
    this->ThresholdFunction = &vtkThresholdPolyData::Upper;
    this->ComponentMode = VTK_COMPONENT_MODE_USE_SELECTED;
    this->SelectedComponent = 0;
    this->PointsDataType = VTK_FLOAT;

    // by default process active point scalars
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
            vtkDataSetAttributes::SCALARS);
}

vtkThresholdPolyData::~vtkThresholdPolyData() {
}

// Criterion is cells whose scalars are less or equal to lower threshold.

void vtkThresholdPolyData::ThresholdByLower(double lower) {
    if (this->LowerThreshold != lower ||
            this->ThresholdFunction != &vtkThresholdPolyData::Lower) {
        this->LowerThreshold = lower;
        this->ThresholdFunction = &vtkThresholdPolyData::Lower;
        this->Modified();
    }
}

// Criterion is cells whose scalars are greater or equal to upper threshold.

void vtkThresholdPolyData::ThresholdByUpper(double upper) {
    if (this->UpperThreshold != upper ||
            this->ThresholdFunction != &vtkThresholdPolyData::Upper) {
        this->UpperThreshold = upper;
        this->ThresholdFunction = &vtkThresholdPolyData::Upper;
        this->Modified();
    }
}

// Criterion is cells whose scalars are between lower and upper thresholds.

void vtkThresholdPolyData::ThresholdBetween(double lower, double upper) {
    if (this->LowerThreshold != lower || this->UpperThreshold != upper ||
            this->ThresholdFunction != &vtkThresholdPolyData::Between) {
        this->LowerThreshold = lower;
        this->UpperThreshold = upper;
        this->ThresholdFunction = &vtkThresholdPolyData::Between;
        this->Modified();
    }
}

int vtkThresholdPolyData::RequestData(
        vtkInformation *vtkNotUsed(request),
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector) {
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and ouptut
    vtkDataSet *input = vtkDataSet::SafeDownCast(
            inInfo->Get(vtkDataObject::DATA_OBJECT()));
    //vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    //  outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkPolyData *output = vtkPolyData::SafeDownCast(
            outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType cellId, newCellId;
    vtkIdList *cellPts, *pointMap;
    vtkIdList *newCellPts;
    vtkCell *cell;
    vtkPoints *newPoints;
    int i, ptId, newId, numPts;
    int numCellPts;
    double x[3];
    vtkPointData *pd = input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *cd = input->GetCellData(), *outCD = output->GetCellData();
    int keepCell, usePointScalars;

    vtkDebugMacro( << "Executing threshold filter");

    if (this->AttributeMode != -1) {
        vtkErrorMacro( << "You have set the attribute mode on vtkThresholdPolyData. This method is deprecated, please use SetInputArrayToProcess instead.");
        return 1;
    }

    vtkDataArray *inScalars = this->GetInputArrayToProcess(0, inputVector);

    if (!inScalars) {
        vtkDebugMacro( << "No scalar data to threshold");
        return 1;
    }

    outPD->CopyGlobalIdsOn();
    outPD->CopyAllocate(pd);
    outCD->CopyGlobalIdsOn();
    outCD->CopyAllocate(cd);

    numPts = input->GetNumberOfPoints();
    output->Allocate(input->GetNumberOfCells());

    newPoints = vtkPoints::New();
    newPoints->SetDataType(this->PointsDataType);
    newPoints->Allocate(numPts);

    pointMap = vtkIdList::New(); //maps old point ids into new
    pointMap->SetNumberOfIds(numPts);
    for (i = 0; i < numPts; i++) {
        pointMap->SetId(i, -1);
    }

    newCellPts = vtkIdList::New();

    // are we using pointScalars?
    usePointScalars = (inScalars->GetNumberOfTuples() == numPts);

    // Check that the scalars of each cell satisfy the threshold criterion
    for (cellId = 0; cellId < input->GetNumberOfCells(); cellId++) {
        cell = input->GetCell(cellId);
        cellPts = cell->GetPointIds();
        numCellPts = cell->GetNumberOfPoints();

        if (usePointScalars) {
            if (this->AllScalars) {
                keepCell = 1;
                for (i = 0; keepCell && (i < numCellPts); i++) {
                    ptId = cellPts->GetId(i);
                    keepCell = this->EvaluateComponents(inScalars, ptId);
                }
            } else {
                keepCell = 0;
                for (i = 0; (!keepCell) && (i < numCellPts); i++) {
                    ptId = cellPts->GetId(i);
                    keepCell = this->EvaluateComponents(inScalars, ptId);
                }
            }
        } else //use cell scalars
        {
            keepCell = this->EvaluateComponents(inScalars, cellId);
        }

        if (numCellPts > 0 && keepCell) {
            // satisfied thresholding (also non-empty cell, i.e. not VTK_EMPTY_CELL)
            for (i = 0; i < numCellPts; i++) {
                ptId = cellPts->GetId(i);
                if ((newId = pointMap->GetId(ptId)) < 0) {
                    input->GetPoint(ptId, x);
                    newId = newPoints->InsertNextPoint(x);
                    pointMap->SetId(ptId, newId);
                    outPD->CopyData(pd, ptId, newId);
                }
                newCellPts->InsertId(i, newId);
            }
            newCellId = output->InsertNextCell(cell->GetCellType(), newCellPts);
            outCD->CopyData(cd, cellId, newCellId);
            newCellPts->Reset();
        } // satisfied thresholding
    } // for all cells

    vtkDebugMacro( << "Extracted " << output->GetNumberOfCells()
            << " number of cells.");

    // now clean up / update ourselves



    pointMap->Delete();
    newCellPts->Delete();

    output->SetPoints(newPoints);
    newPoints->Delete();

    output->Squeeze();
    return 1;
}

int vtkThresholdPolyData::EvaluateComponents(vtkDataArray *scalars, vtkIdType id) {
    int keepCell = 0;
    int numComp = scalars->GetNumberOfComponents();
    int c;

    switch (this->ComponentMode) {
        case VTK_COMPONENT_MODE_USE_SELECTED:
            c = (this->SelectedComponent < numComp) ? (this->SelectedComponent) : (0);
            keepCell =
                    (this->*(this->ThresholdFunction))(scalars->GetComponent(id, c));
            break;
        case VTK_COMPONENT_MODE_USE_ANY:
            keepCell = 0;
            for (c = 0; (!keepCell) && (c < numComp); c++) {
                keepCell =
                        (this->*(this->ThresholdFunction))(scalars->GetComponent(id, c));
            }
            break;
        case VTK_COMPONENT_MODE_USE_ALL:
            keepCell = 1;
            for (c = 0; keepCell && (c < numComp); c++) {
                keepCell =
                        (this->*(this->ThresholdFunction))(scalars->GetComponent(id, c));
            }
            break;
    }
    return keepCell;
}

// Return the method for manipulating scalar data as a string.

const char *vtkThresholdPolyData::GetAttributeModeAsString(void) {
    if (this->AttributeMode == VTK_ATTRIBUTE_MODE_DEFAULT) {
        return "Default";
    } else if (this->AttributeMode == VTK_ATTRIBUTE_MODE_USE_POINT_DATA) {
        return "UsePointData";
    } else {
        return "UseCellData";
    }
}

// Return a string representation of the component mode

const char *vtkThresholdPolyData::GetComponentModeAsString(void) {
    if (this->ComponentMode == VTK_COMPONENT_MODE_USE_SELECTED) {
        return "UseSelected";
    } else if (this->ComponentMode == VTK_COMPONENT_MODE_USE_ANY) {
        return "UseAny";
    } else {
        return "UseAll";
    }
}

int vtkThresholdPolyData::FillInputPortInformation(int, vtkInformation *info) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
}

void vtkThresholdPolyData::PrintSelf(ostream& os, vtkIndent indent) {
    this->Superclass::PrintSelf(os, indent);

    os << indent << "Attribute Mode: " << this->GetAttributeModeAsString() << endl;
    os << indent << "Component Mode: " << this->GetComponentModeAsString() << endl;
    os << indent << "Selected Component: " << this->SelectedComponent << endl;

    os << indent << "All Scalars: " << this->AllScalars << "\n";
    if (this->ThresholdFunction == &vtkThresholdPolyData::Upper) {
        os << indent << "Threshold By Upper\n";
    }
    else if (this->ThresholdFunction == &vtkThresholdPolyData::Lower) {
        os << indent << "Threshold By Lower\n";
    }
    else if (this->ThresholdFunction == &vtkThresholdPolyData::Between) {
        os << indent << "Threshold Between\n";
    }

    os << indent << "Lower Threshold: " << this->LowerThreshold << "\n";
    os << indent << "Upper Threshold: " << this->UpperThreshold << "\n";
    os << indent << "DataType of the output points: "
            << this->PointsDataType << "\n";
}
