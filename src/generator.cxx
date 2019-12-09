#include <vtkTessellatedBoxSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

int main(){
	vtkNew<vtkTessellatedBoxSource> src;
	src->SetLevel(30);
	src->SetBounds(0,200,-20,20,-20,20);
	//src->SetLevel(1);
	//src->SetBounds(-20,5,-30,30,-30,30);
	vtkNew<vtkPolyDataWriter> wri;
	wri->SetFileName("out.vtk");
	wri->SetInputConnection(src->GetOutputPort());
	wri->Write();
	return EXIT_SUCCESS;
}

