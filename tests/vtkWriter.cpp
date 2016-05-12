#include "vtkWriter.h"
#include "node.h"


#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDirectory.h>
#include <vtkLine.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>

VtkWriter::VtkWriter(const char* file_name)
{
    this->file_name = file_name;
}

int VtkWriter::writeCase(paraFEM::FemCasePtr c, double coordSysSize)
{
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();  // data container
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();  // nodes
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();  // truss elements
    vtkSmartPointer<vtkFloatArray> cellData = vtkSmartPointer<vtkFloatArray>::New();  // truss forces
    vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();  // polygons elements (t3/t4)

    //set nodes
    for (auto node: c->nodes)
    {
        paraFEM::Vector3 pos = node->position;
        points->InsertNextPoint(pos.x(), pos.y(), pos.z());
    }

    //set elements + forces
    for (auto element: c->elements)
    {
        vector<int> numbers = element->getNr();
        if (numbers.size() == 2) { // truss element
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            line->GetPointIds()->SetId(0, numbers[0]);
            line->GetPointIds()->SetId(1, numbers[1]);
            lines->InsertNextCell(line);
            paraFEM::TrussPtr truss = std::dynamic_pointer_cast<paraFEM::Truss>(element);
            cellData->InsertNextValue(truss->stress);  //scalar
        } else if(numbers.size() > 2) { // membrane element (3/ 4 nodes)
            paraFEM::MembranePtr membrane = std::dynamic_pointer_cast<paraFEM::Membrane>(element);
            if (membrane and coordSysSize <=0)
            {
                 cellData->InsertNextValue(pow(membrane->stress.x(),2) + pow(membrane->stress.y(),2));
            }
            vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
            for (int nr: numbers) {
                polygon->GetPointIds()->InsertNextId(nr);
            }
            polygons->InsertNextCell(polygon);
        }

    }
    if (coordSysSize > 0)
    {
        int count = c->nodes.size();
        paraFEM::Vector3 center; //count
        paraFEM::Vector3 x;      // count +1
        paraFEM::Vector3 y;      // count +2
        paraFEM::Vector3 z;      // count +3
        for (auto element: c->elements)
        {
            paraFEM::MembranePtr membrane = std::dynamic_pointer_cast<paraFEM::Membrane>(element);
            if (membrane)
            {
                center = membrane->center;
                x = center + membrane->coordSys.t1 * coordSysSize;
                y = center + membrane->coordSys.t2 * coordSysSize;
                z = center + membrane->coordSys.n * coordSysSize;
                for (auto pos: vector<paraFEM::Vector3>{center, x, y, z})
                    points->InsertNextPoint(pos.x(), pos.y(), pos.z());
                for (int i = 0; i < 3; i++)
                {
                    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                    line->GetPointIds()->SetId(0, count);
                    line->GetPointIds()->SetId(1, count + i + 1);
                    lines->InsertNextCell(line);
                }
                count += 4;
            }
        }
    }

    polydata->SetPoints(points);
    polydata->SetLines(lines);
    polydata->SetPolys(polygons);
    polydata->GetCellData()->SetScalars(cellData);

    writer->SetFileName((file_name + to_string(count) + string(".vtk")).c_str());
    writer->SetInputData(polydata);

    // create directory if needed
    size_t found = file_name.find_last_of("/\\");
    if (found != -1){
        string dirname = file_name.substr(0, found);
        vtkDirectory *directory = vtkDirectory::New();
        directory->MakeDirectory(dirname.c_str());
    }

    int done  = writer->Write();
    count ++;

    return done;
}