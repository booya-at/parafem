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

namespace paraFEM{

VtkWriter::VtkWriter(const char* file_name)
{
    this->file_name = file_name;
}

int VtkWriter::writeCase(paraFEM::FemCasePtr c, double coordSysSize)
{
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();  // data container
    vtkSmartPointer<vtkPolyData> polydata_2 = vtkSmartPointer<vtkPolyData>::New();  // data container
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();  // nodes
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();  // truss elements
    vtkSmartPointer<vtkFloatArray> cellData = vtkSmartPointer<vtkFloatArray>::New();  // truss forces
    vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();  // polygons elements (t3/t4)
    
    cellData->SetNumberOfComponents(3);
    cellData->SetName("stress");
    
    std::vector<paraFEM::TrussPtr> truss_elements;
    std::vector<paraFEM::MembranePtr> membrane_elements;
    
    for (auto element: c->elements)
    {
        paraFEM::TrussPtr t = std::dynamic_pointer_cast<paraFEM::Truss>(element);
        paraFEM::MembranePtr m = std::dynamic_pointer_cast<paraFEM::Membrane>(element);
        if (t)
            truss_elements.push_back(t);
        if (m)
            membrane_elements.push_back(m);
    }

    //set nodes
    int counter = 0;
    if (coordSysSize > 0)
    {
        paraFEM::Vector3 center; //counter
        paraFEM::Vector3 x;      // counter +1
        paraFEM::Vector3 y;      // counter +2
        paraFEM::Vector3 z;      // counter +3
        for (auto element: membrane_elements)
        {
            center = element->center;
            x = center + element->coordSys.t1 * coordSysSize;
            y = center + element->coordSys.t2 * coordSysSize;
            z = center + element->coordSys.n * coordSysSize;
            for (auto pos: vector<paraFEM::Vector3>{center, x, y, z})
                points->InsertNextPoint(pos.x(), pos.y(), pos.z());
            for (int i = 0; i < 3; i++)
            {
                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                line->GetPointIds()->SetId(0, counter);
                line->GetPointIds()->SetId(1, counter + i + 1);
                lines->InsertNextCell(line);
                
                cellData->InsertNextTuple3(0, 0, 0);
            }
            counter += 4;
        }
    }
    
    
    for (auto node: c->nodes)
    {
        paraFEM::Vector3 pos = node->position;
        points->InsertNextPoint(pos.x(), pos.y(), pos.z());
    }

    for (auto element: truss_elements)
    {
        vector<int> numbers = element->getNr();
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, numbers[0] + counter);
        line->GetPointIds()->SetId(1, numbers[1] + counter);
        lines->InsertNextCell(line);
        paraFEM::Vector3 stress = element->getStress();
        cellData->InsertNextTuple3(stress.x(), stress.y(), stress.z());  //scalar
    }

    for (auto element: membrane_elements)
    {
        vector<int> numbers = element->getNr();
        vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
        for (int nr: numbers) {
            polygon->GetPointIds()->InsertNextId(nr + counter);
        }
        polygons->InsertNextCell(polygon);
        paraFEM::Vector3 stress = element->getStress();
        cellData->InsertNextTuple3(stress.x(), stress.y(), stress.z());
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

}