/* ###########################################################
   # PBD algorithm implementation                            #
   # Thesis: Muscle Fibres Deformation using Particle System #
   # Author: Martin Cervenka                                 #
   # Version: 5/2019                                         #
   ###########################################################

	This module integrates VTK with PBD algorithm and
	provides basic functionality as a demo application
*/
#include <vtkCylinderSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPolyDataWriter.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkSmartPointer.h>
#include <vtkCubeSource.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProgrammableFilter.h>
#include <vtkCallbackCommand.h>
#include <vtkAppendPolyData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCleanPolyData.h>
#include <vtkCornerAnnotation.h>
#include <vtkTextProperty.h>
#include <vtkSliderWidget.h>
#include <vtkProperty2D.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <memory>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <iostream>

#include "pbd.h"
/* Typedef for better readability */
typedef vtkSliderRepresentation2D vtkSlider;
/* Actual iteration */
static int iteration = 0;
/* Width of the window */
constexpr size_t WIDTH = 600;
/* Height of the window */
constexpr size_t HEIGHT = 600;
/* Parameters for one algorithm iteration */
struct vtkParams{
	/* pbd algorithm pointer */
	pbd* algorithm;
	/* data for the algorithm */
	data* dta;
	/* VTK deformed polydata */
	vtkSmartPointer<vtkPolyData> polydata;
	/* VTK Fiter with output polydata */
	vtkSmartPointer<vtkProgrammableFilter> filter;
	/* Label with FPS and iteration number */
	vtkSmartPointer<vtkCornerAnnotation> fps_string;
	/* Vector with frame numbers, in which we get screens */
	std::vector<size_t> export_states;
	/* and current window */
	vtkSmartPointer<vtkRenderWindow> window;
};
/*
	Function is a callback for every slider change in the window
	=> caller Calling object (slider)
	=> client_data Associated data to the slider (double pointer)
*/
void SliderCallback(vtkObject* caller,long unsigned,void* client_data,void*){
	/* We know, that the caller is a slider widget */
	vtkSliderWidget* wid = reinterpret_cast<vtkSliderWidget*>(caller);
	/* Then we need the slider object from the slider widget */
	vtkSlider* slid = reinterpret_cast<vtkSlider*>(wid->GetRepresentation());
	/* We retrieve value from the slider */
	double value = slid->GetValue();
	/* And finally set the client data to the given value */
	*reinterpret_cast<double*>(client_data) = value;
}
/*
	Function makes screenshot of the window and saves it into file
	=> pth Output file path
	=> renderWindow Window to shot
*/
void screenshot(const char* pth, vtkSmartPointer<vtkRenderWindow> renderWindow){
	/* first we create filter transforming window to image */
	vtkNew<vtkWindowToImageFilter> windowToImageFilter;
	/* we declare which window we want */
	windowToImageFilter->SetInput(renderWindow);
	/* then, we specify image type */
	windowToImageFilter->SetInputBufferTypeToRGBA();
	/* we do not want to read front buffer */
	windowToImageFilter->ReadFrontBufferOff();
	/* now, we can update this filter */
	windowToImageFilter->Update();
	/* then, we prepare filter which saves image to the output file */
	vtkNew<vtkPNGWriter> writer;
	/* so we specify the path */
	writer->SetFileName(pth);
	/* and source data connection */
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	/* finally, we create the file */
	writer->Write();
}
/*
	This method is called at regular intervals - for each algorithm iteration
	=> caller Calling object (Interactor)
	=> client_data Data associated with the event (one iteration filter)
*/
void TimerCallback(vtkObject* caller,long unsigned,void* client_data,void*){
	/* Retreives the one iteration filter */
	vtkSmartPointer<vtkProgrammableFilter> pf =
		static_cast<vtkProgrammableFilter*>(client_data);
	/* And interactor */
	vtkRenderWindowInteractor* iren =
		static_cast<vtkRenderWindowInteractor*>(caller);
	/* We have to tell filter it has to be recalculated
		(it has to make new iteration) */
	pf->Modified();
	/* And finally we render the change */
	iren->Render();
}
/*
	Function which frees memory at the end
	=> data Data for deallocation (struct vtkParams)
*/
void DeleteArg(void* data){
	/* At first we get the vtkParams struct from void pointer */
	struct vtkParams* param = static_cast<vtkParams*>(data);
	/* And then we free the data structure */
	delete param->dta;
	/* PBD algorithm */
	delete param->algorithm;
	/* And finally the whole param, too */
	delete param;
}
/*
	One iteration of the algorithm
	=> data Iteration parameters (struct vtkParams)
*/
void OneIteration(void* data){
	/* We get the parameters from the void pointer (nothing fancy, just cast) */
	struct vtkParams* param = static_cast<vtkParams*>(data);
	/* Then we need to get actual time (for FPS and time measurement) */
	auto t1 = std::chrono::high_resolution_clock::now();
	/* We try to do one PBD algrithm iteration */
	if(!param->algorithm->execute(iteration)){
		/* it it was successful, we can increment the actual iteration value */
		iteration++;
	} 
	/* We get actual time (for time measurement) */
	CAPTURE_TIME(t2);
	/* For the output poly data */
	vtkPolyData* out = param->filter->GetPolyDataOutput();
	/* We just copy the input data, which has been changed by PBD algorithm */
	out->DeepCopy(param->polydata);
	/* Now we get the time after PBD */
	auto t3 = std::chrono::high_resolution_clock::now();
	/* And calculate spend time */
	auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>
		(t3-t1).count();
	/* Then, we print some information for time measurement */
	PRINT_DIFF("Coping data:              ",cpy,t3,t2);
	PRINT_DIFF("Total time (approx.):     ",tmp,t3,t1);
	/* If we want to lock FPS at some level */
#ifdef FPS_LIMIT
	/* We calculate number of nanoseconds from FPS */
	constexpr int lock_fps = 1e9/FPS_LIMIT;
	/* If algoritm was faster than the limit, we have to wait certain time */
	if(diff<lock_fps){
		/* We prepare structure for sleep time information */
		struct timespec sleep_time,slp;
		/* Then, we set sleep time in second to zero
			(assuming FPS won't be <1) */
		sleep_time.tv_sec = 0;
		/* Next, we set nanoseconds to fill the remaining time */
		sleep_time.tv_nsec = (lock_fps-diff); 
		/* We take a break */
		nanosleep(&sleep_time,&slp); 
		/* And then we have to rewrite the ending time */
		auto t4 = std::chrono::high_resolution_clock::now();
		/* The same with time difference */
		diff = std::chrono::duration_cast<std::chrono::nanoseconds>
			(t4-t1).count();
	}
#endif
	/* And finally the label with FPS and actual iteration is set */
	param->fps_string->SetText(vtkCornerAnnotation::UpperLeft,
		(std::string("ITER: ")+std::to_string(iteration)+
		std::string(", FPS: ")+
		std::to_string(static_cast<double>(1e9/diff))).c_str()
	);
	/* Last what we can is to create screenshot. We get number of frames */
	size_t len = param->export_states.size();
	/* and prepare buffer for output filename */
	char buff[50];
	/* for each screenshot */
	for(size_t a=0; a<len; a++){
		/* we get frame, when it should be screened */
		size_t val = param->export_states[a];
		/* it frame ID is the same as current interation */
		if(static_cast<size_t>(iteration)==val){
			/* we prepare filename */
			snprintf(buff,50,"screen_%04lu.png",val);
			/* and do the screenshot */
			screenshot(buff,param->window);
		}
	}
}

/*
	This function creates one slider
	=> min Minimal slider value
	=> max Maximal slider value
	=> value Actual value/value which will be set with this slider
	=> title Title of the slider
	=> x1 X coordinate of the beginning of the slider
	=> y1 Y coordinate of the beginning of the slider
	=> x2 X coordinate of the end of the slider
	=> y2 Y coordinate of the end of the slider
	=> interactor VTK interactor, with which we can control the slider
*/
vtkSmartPointer<vtkSliderWidget> get_slider(
	double min,double max,double* value,const char* title,size_t x1,size_t y1,
	size_t x2, size_t y2,
	vtkSmartPointer<vtkRenderWindowInteractor> interactor
){
	/* At first, we create the slider */
	vtkNew<vtkSlider> slider;
	/* Then, we set minimum */
	slider->SetMinimumValue(min);
	/* and maximum value */
	slider->SetMaximumValue(max);
	/* We set the value */
	slider->SetValue(*value);
	/* And title */
	slider->SetTitleText(title);
	/* We switch the coordinate system to dislay for
		the beginning point of the slider */
	slider->GetPoint1Coordinate()->SetCoordinateSystemToDisplay();
	/* And set the coordinates according to the function input */
	slider->GetPoint1Coordinate()->SetValue(x1,y1);
	/* And the same for the second point, we switch the coordinate system to
		dislay for the ending point of the slider */
	slider->GetPoint2Coordinate()->SetCoordinateSystemToDisplay();
	/* And set the coordinates according to the function input */
	slider->GetPoint2Coordinate()->SetValue(x2,y2);
	/* The slider is red */
	slider->GetSliderProperty()->SetColor(1,0,0);
	/* Title is dark green */
	slider->GetTitleProperty()->SetColor(0,.8,0);
	/* Label is dark green too */
	slider->GetLabelProperty()->SetColor(0,.8,0);
	/* The line after which the slider slides is light gray */
	slider->GetTubeProperty()->SetColor(.7,.7,.7);
	/* And the caps have the same color (light gray) */
	slider->GetCapProperty()->SetColor(.7,.7,.7);
	/* Then, we turn off the shadow for the title (looks bad on screenshots) */
	slider->GetTitleProperty()->ShadowOff();
	/* and turn off the shadow for the label, too (same reason) */
	slider->GetLabelProperty()->ShadowOff();
	/* Original line width is too thick, so we have to adjust it */
	slider->GetTubeProperty()->SetLineWidth(0.1);
	/* same with caps */
	slider->GetCapProperty()->SetLineWidth(0.1);
	/* Next step is create the widget for the slider */
	vtkNew<vtkSliderWidget> widget;
	/* We set the interactor (else we won't be able to control the widget) */
	widget->SetInteractor(interactor);
	/* Then, we assign the slider to the widget */
	widget->SetRepresentation(slider);
	/* We allow animation for the slider widget */
	widget->SetAnimationModeToAnimate();
	/* And it has to be enabled */
	widget->EnabledOn();
	/* Finally it is necessary to assign a callback to be useful */
	vtkNew<vtkCallbackCommand> sliderCallback;
	/* We pass the function pointer of the callback */
	sliderCallback->SetCallback(SliderCallback);
	/* set the client data (pointer to the value, which has to be changed) */
	sliderCallback->SetClientData(value);
	/* and the operation has to be bound with the widget */
	widget->AddObserver(vtkCommand::InteractionEvent,sliderCallback);
	/* That's all, we can return the widget */
	return widget;
}
/*
	Function creates all sliders and puts them to the window
	=> pbd_algorithm PBD algorithm - it is required because there
		are variables, which must respond to the slider change
	=> renderWindowInteractor Interactor, we have to be able
		to interact with them
	=> slider0 Created slider - "Velocity damp"
	=> slider1 Created slider - "Bending"
	=> slider2 Created slider - "Iteration / frame"
	=> slider3 Created slider - "Pressure"
	=> slider4 Created slider - "Anisotropy"
	=> slider5 Created slider - "Gravity"
	=> slider6 Created slider - "Pause 0/1"
*/
void create_sliders(
	pbd* pbd_algorithm,
	vtkRenderWindowInteractor* renderWindowInteractor,
	vtkSmartPointer<vtkSliderWidget>& slider0,
	vtkSmartPointer<vtkSliderWidget>& slider1,
	vtkSmartPointer<vtkSliderWidget>& slider2,
	vtkSmartPointer<vtkSliderWidget>& slider3,
	vtkSmartPointer<vtkSliderWidget>& slider4,
	vtkSmartPointer<vtkSliderWidget>& slider5,
	vtkSmartPointer<vtkSliderWidget>& slider6
){
	/* First slider sets value of "Velocity damp", value
		is in interval <0;1>, and position is [20;100], [200;100] */
	slider0 = get_slider(
		0,1,&pbd_algorithm->m_damp,"Velocity damp",
		20,100,200,100,renderWindowInteractor
	);
	/* Second slider sets bending factor */
	slider1 = get_slider(
		0,1,&pbd_algorithm->m_bending,"Bending",
		20,200,200,200,renderWindowInteractor
	);
	/* Third slider sets number of iterations per one frame */
	slider2 = get_slider(
		1,100,&pbd_algorithm->m_constr,"Iteration / frame",
		20,300,200,300,renderWindowInteractor
	);
	/* Fourth slider is for pressure, it changes volume
		of the muscle */
	slider3 = get_slider(
		0.5,2,&pbd_algorithm->m_pressure,"Pressure",
		20,400,200,400,renderWindowInteractor
	);
	/* Next slider sets anisotropy */
	slider4 = get_slider(
		0,2,&pbd_algorithm->m_anisotropy,"Anisotropy",
		20,500,200,500,renderWindowInteractor
	);
	/* The last but one slider sets gravity */
	slider5 = get_slider(
		0,5,&pbd_algorithm->m_gravity,"Gravity",
		280,100,400,100,renderWindowInteractor
	);
	/* And last is for pause (this is tricky one, because
		it is only binary, it is pause if the value is
		below 0.5) */
	slider6 = get_slider(
		0,1,&pbd_algorithm->m_pause,"Pause 0/1",
		480,100,600,100,renderWindowInteractor
	);
}
/*
	This function parses input parameters
	=> argc Parameter count
	=> argv Parameter values
	<=> len Point count for each model (single muscle/bone)
	<=> classes Model count in each category (muscles/bones)
	<=> scenario ID of the dynamic scenario
	<=> fibres Number of fibre datasets in the simulation
	<=> fibres_append Filter which connects fibre datasets together
	<= All models connected together
*/
vtkSmartPointer<vtkAppendPolyData> parse_input(
	int argc,
	char** argv,
	std::vector<size_t>& len,
	std::vector<size_t>& classes,
	std::vector<size_t>& export_states,
	std::vector<double>& parameters,
	size_t& scenario,
	size_t& fibres,
	vtkSmartPointer<vtkAppendPolyData> fibres_append
){
	/* Model count in current class (bones/muscles etc.) */
	size_t count = 0;
	/* Do we set scenario with next argument or not? */
	bool scenario_state = false;
	/* Structure for all data */
	vtkNew<vtkAppendPolyData> append;
	/* Next, we specify if we are in export argument or not */
	bool export_state = false;
	bool parameters_state = false;
	/* At first, we iterate over all arguments */
	for(int a=1; a<argc; a++){
		/* If the argument defines muscle class */
		if(!strcmp(argv[a],"--muscles")){
			/* Then we just ignore it (there is no point writing
				this argument) */
			continue;
		}
		/* If the argument is "bones" */
		if(!strcmp(argv[a],"--bones")){
			/* We save number of muscle models into output vector */
			classes.push_back(count);
			/* Then reset the model counter (now it will count bones) */
			count=0;
			/* And continue with next iteration */
			continue;
		}
		/* If fibres follows */
		if(!strcmp(argv[a],"--fibres")){
			/* Then we just switch the variable, which indicates fibres */
			fibres=1;
			/* And continue with next argument */
			continue;
		}
		/* If we want to specify scenario ID */
		if(!strcmp(argv[a],"--scenario")){
			/* We switch the scenario flag */
			scenario_state=true;
			/* And skip this argument */
			continue;
		}
		/* If this parameter declares "export" section */
		if(!strcmp(argv[a],"--export")){
			/* then we switch to export state */
			export_state=true;
			/* and nothing else is required */
			continue;
		}
		if(!strcmp(argv[a],"--parameters")){
			parameters_state=true;
			continue;
		}
		if(parameters_state){
			parameters.push_back(atof(argv[a]));
		/* If we are selecting the scenario right now */
		}else if(scenario_state){
			/* Next argument won't select scenario */
			scenario_state=false;
			/* And this scenario is parsed */
			scenario = atoi(argv[a]);
		/* In case we want export screenshots */
		}else if(export_state){
			/* we save iteration ID */
			export_states.push_back(atoi(argv[a]));
		/* In case we want fibres */
		}else if(fibres){
			/* We prepare the VTK reader for the fibre dataset */
			vtkNew<vtkGenericDataObjectReader> fibres_data;
			/* Then we specify the path to the dataset file */
			fibres_data->SetFileName(argv[a]);
			/* And put the filter into the filter chain */
			fibres_append->AddInputConnection(fibres_data->GetOutputPort());
			/* We have loaded one more dataset, so we increment the variable */
			fibres++;
		/* If muscle or bone model is being loaded */
		}else{
			/* We declare VTK reader for single model */	
			vtkNew<vtkGenericDataObjectReader> muscle_data;
			/* We also prepare the color array */
			vtkNew<vtkUnsignedCharArray> colors;
			/* which will have 3 components (RGB) */
			colors->SetNumberOfComponents(3);
			/* and it's name will be "Colors" */
			colors->SetName("Colors");
			/* We set the path to the model */
			muscle_data->SetFileName(argv[a]);
			/* get the polydata from the reader */
			vtkSmartPointer<vtkPolyData> data=muscle_data->GetPolyDataOutput();
			/* and update the data (so we load the model) */
			muscle_data->Update();
			/* Color initialization follows, at first, we get the loop count
				(number of point, which will have certain color) */
			size_t points = data->GetNumberOfPoints();
			/* In the meantime we save point number into the vector */
			len.push_back(points);
			/* For each point in the curent model */
			for(size_t a=0; a<points; a++){
				/* For each color component */
				for(size_t b=0; b<3; b++){
					/* set its value to 200 (it becomes light gray) */
					colors->InsertNextValue(100);
				}
			}
			/* Then we assign created color to the model */
			data->GetPointData()->SetScalars(colors);
			/* add current data to the rest */
			append->AddInputConnection(muscle_data->GetOutputPort());
			/* and increment number of models in current class */
			count++;
		}
	}
	/* Finally, we save last model count in the class count vector */
	classes.push_back(count);
	/* Then we update filter, which merges all loaded data */
	append->Update();
	/* and return it */
	return append;
}
/*
	Main function
	=> argc Number of arguments
	=> argv Argument values
*/
int main(int argc,char** argv){
	/* At first, we create the PBD algorithm (delete in DeleteArg function) */
	pbd* pbd_algorithm = new pbd();
	/* Then, we create window */
	vtkNew<vtkRenderWindow> renderWindow;
	/* renderer */
	vtkNew<vtkRenderer> renderer;
	/* one iteration VTK algorithm (which just calls PBD) */
	vtkNew<vtkProgrammableFilter> algo;
	/* mapper for the deforming models */
	vtkNew<vtkPolyDataMapper> mapper;
	/* and actor */
	vtkNew<vtkActor> actor;
	/* Next, we need interactor for sliders */
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	/* timer callback for automatic PBD iterations */
	vtkNew<vtkCallbackCommand> timerCallback;
	/* top left corner text (FPS and iteration) */
	vtkNew<vtkCornerAnnotation> textActor;
	/* and filter which merges all fibre direction data */
	vtkNew<vtkAppendPolyData> fibres_append;
	/* We prepare "empty" top left corner text */
	textActor->SetText(vtkCornerAnnotation::UpperLeft,"FPS: -");
	/* We declare its size to 10 (to be well visible) */
	textActor->GetTextProperty()->SetFontSize(10);
	/* and set its color to green (visible on black background) */
	textActor->GetTextProperty()->SetColor(0,.8,0);
	/* Next step is connection render to window */
	renderWindow->AddRenderer(renderer);
	/* and iteractor to window too */
	renderWindowInteractor->SetRenderWindow(renderWindow);
/* If we want to see HUD (sliders and top left text) */
	/* We need vector, which contains number of points in each loaded model */
	std::vector<size_t> len;
	/* and a vector which contains number of models in each class
		(muscles/bones etc.) */
	std::vector<size_t> classes;
	/* and then we save all frame indices, where screenshot should occur */
	std::vector<size_t> export_states;
	std::vector<double> parameters;
	/* Next, we prepare variables for scenario ID and
		number of fibres datasets */
	size_t scenario=0,fibres=0;
	/* At this point, we load all input data accoring to the parameters from
		the command line */
	vtkSmartPointer<vtkAppendPolyData> append = parse_input(
		argc,argv,len,classes,export_states,parameters,scenario,fibres,fibres_append
	);
	pbd_algorithm->set_parameters(parameters);
#ifndef NO_HUD
	/* In this case we add annotation */
	renderer->AddActor2D(textActor);
	/* Prepare smart pointer for every slider */
	vtkSmartPointer<vtkSliderWidget> s0,s1,s2,s3,s4,s5,s6;
	/* and create them */
	create_sliders(pbd_algorithm,renderWindowInteractor,s0,s1,s2,s3,s4,s5,s6);
#endif
	/* We need array, which will contain fibres position in doubles */
	vtkNew<vtkDoubleArray> double_fibres;
	/* If fibres has been loaded */
	if(fibres>1){
		/* We update the fibres filter */
		fibres_append->Update();
		/* and copy positions to the double array (we dont know apriori, if the
			input data are in FLOATs or DOUBLEs, that's why this is done) */
		double_fibres->DeepCopy(
			fibres_append->GetOutput()->GetPoints()->GetData()
		);
	}
	/* Next step is to create a new double array for all model points
		coordinates */
	vtkNew<vtkDoubleArray> double_arr;
	/* We copy the original coordinates (don't know if it was FLOATs or DOUBLEs)
		into the new DOUBLE array, so we know the data type now */
	double_arr->DeepCopy(append->GetOutput()->GetPoints()->GetData());
	/* Then we need to create control structure for the PBD algorithm */
	vtkParams* param = new vtkParams();
	/* which is then filled with all data */
	param->polydata = append->GetOutput();
	/* but coordinates are changed to DOUBLEs */
	param->polydata->GetPoints()->SetData(double_arr);
	/* We pass the reference to the top left label, too */
	param->fps_string = textActor;
	/* pass export frame IDs too */
	param->export_states = export_states;
	/* pass window */
	param->window = renderWindow;
	/* and create data structure (independent of VTK library) */
	param->dta = new data (
		reinterpret_cast<double*>(double_arr->GetVoidPointer(0)),
		param->polydata->GetPoints()->GetNumberOfPoints(),
		param->polydata->GetPolys()->GetPointer(),
		param->polydata->GetNumberOfPolys(),
		reinterpret_cast<unsigned char*>(
			param->polydata->GetPointData()->GetScalars()->GetVoidPointer(0)
		),
		param->polydata->GetPointData()->GetScalars()->GetNumberOfComponents(),
		reinterpret_cast<double*>(double_fibres->GetVoidPointer(0)),
		(fibres>1)?
			fibres_append->GetOutput()->GetPoints()->GetNumberOfPoints():0,
		(fibres>1)?
			fibres_append->GetOutput()->GetLines()->GetPointer():nullptr,
		(fibres>1)?
			fibres_append->GetOutput()->GetNumberOfLines():0
	);
	/* Then, we init PBD algorithm */
	pbd_algorithm->init(
		param->dta->create_graph(len,classes,1000),param->dta,scenario
	);
	/* pass it into the parameter structure */
	param->algorithm = pbd_algorithm;
	/* and pass the one iteration filter, too */
	param->filter = algo;
	/* We then initialize one iteration filter, which will be executed
		by timer */
	algo->SetInputConnection(append->GetOutputPort());
	/* One iteration will call "OneIteration" function */
	algo->SetExecuteMethod(OneIteration,param);
	/* and cleans memory when done */
	algo->SetExecuteMethodArgDelete(DeleteArg);
	/* We have to connect the algorithm into the VTK filter chain */
	mapper->SetInputConnection(algo->GetOutputPort());
	/* Actor has to know it's mapper */
	actor->SetMapper(mapper);
	/* and we set the ambient lighting to 0.1 (it looks smoother) */
	actor->GetProperty()->SetAmbient(0.1);
	/* We assign actor to the renderer */
	renderer->AddActor(actor);
	/* and it's background will be black (let's take care of our eyes) */
	renderer->SetBackground(0,0,0);
	/* Next, we reset the camera position and angle, which calculates
		optimal coordinates to fit everything on screen */
	renderer->ResetCamera();
	/* Window size has to be set */
	renderWindow->SetSize(WIDTH,HEIGHT);
	/* and interactor has to be initialized (now we are able to interact
		with sliders) */
	renderWindowInteractor->Initialize();
	/* Last but not least we create timer callback (to animate simulation) */
	timerCallback->SetCallback(TimerCallback);
	/* Callback needs to know what to do in every iteration */
	timerCallback->SetClientData(algo);
	/* Callback is assigned to the interactor, it will be called now
		whenever interactor needs to */
	renderWindowInteractor->AddObserver(vtkCommand::TimerEvent,timerCallback);
	/* Let's calculate PBD after every 1ms (as fast as possible) */
	renderWindowInteractor->CreateRepeatingTimer(1);
	/* Finally, we show the window */
	renderWindowInteractor->Start();
}
