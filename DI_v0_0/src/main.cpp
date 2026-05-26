// header
#include "_header.hpp"
#include "_main.hpp"
#include <H5Cpp.h>

using namespace std;

int main(int argc, char **argv)
{
	cout << "----- Double Inverse Engagement Program -----" << endl << endl;

	// Avoid HDF5 global atexit finalization crash observed after analysis end.
	try
	{
		H5::H5Library::dontAtExit();
	}
	catch (...)
	{
		// Continue even if the runtime does not support this hook.
	}

	// omp_set_dynamic(1);
	kmp_set_defaults("KMP_BLOCKTIME=0");
	kmp_set_blocktime(0);

	information info;
	int alloc_count = 0, glo_var_count = 0;
	Total_mesh = argc - 3;

	cout << Total_mesh << " input files detected." << endl << endl;

	printf("start Get Option Data\n\n");

	// get option file
	Get_Option(argv[argc - 1], info.opt_files);

	// set option data
	Set_Option_Data(&info);

	// get constant data
	info.c.get_data(info.opt_files[0].c_str(), &info);

	// 処理時間計測
	double time;
	chrono::system_clock::time_point start[2], end[2];
	start[0] = chrono::system_clock::now();
	start[1] = chrono::system_clock::now();

	// check arguments, IGA: input file = 1, S-IGA: input file > 1
	if (Total_mesh <= 0)
		printf("Argument is missing\n\n");
	else if (Total_mesh == 1)
		printf("IGA carried out.(No local mesh)\n\n");
	else if (Total_mesh == 2)
		printf("SS-IGA model data carried out.(%d input files)\n\n", Total_mesh);
	else
	{
		printf("Too many input files. Maximum is 4.\n\n");
		exit(0);
	}

	printf("start Get Input Data\n\n");

	// memory allocation
	printf("before first Allocation\n");
	fflush(stdout);
	Allocation(alloc_count++, &info);

	// Read file 1st time
	for (int i = 0; i < Total_mesh; i++)
	{
		printf("before Get_Input_1 mesh=%d\n", i);
		fflush(stdout);
		Get_Input_1(i, argv[i + 1], &info);
		printf("after Get_Input_1 mesh=%d\n", i);
		fflush(stdout);
	}

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 2nd time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_2(i, argv[i + 1], &info);

	// Read file 3rd time
	printf("\nstart reading vertex data\n");
	Allocation(10, &info);
	Get_Input_3(argv[argc - 2], &info);

	// MAX value
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// INC 等の作成
	printf("\nstart Make INC\n");
	if (Total_mesh >= 2)
		Make_INC_for_local_geometric(&info);
	Make_INC(&info);
	
	// set D_MATRIX_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	printf("start Preprocessing\n");
	if (Total_mesh == 1) // IGA
	{
		Make_Gauss_points(true, &info);
	}
	else if (Total_mesh >= 2) // SS-IGA
	{
		searchOverlappingEle(&info);
		Make_Gauss_points(false, &info);
	}

	// memory allocation
	Allocation(alloc_count++, &info);

	// for output vertex
	if (info.c.CALC_ON_ELE_VERTEX == 1)
	{
		Allocation(9, &info);
		printf("start calc on element vertex\n\n");
		Calc_on_Element_Vertex(&info);

		if (argv[argc - 2] == nullptr)
		{
			printf("vertex data is missing: only vertex output mode\n\n");
			exit(0);
		}
	}

	// check geometry only output
	bool isGeometryOnly = (info.c.GEOMETRY_ONLY_OUTPUT == 1) ? true : false;
	if (isGeometryOnly && info.c.OUTPUT_GLOBAL_PARAMETERS == 1)
	{
		printf("start Make_connectivity\n");
		Allocation(7, &info);
		Allocation(8, &info);
		Make_connectivity(&info);

		printf("Geometry only output mode in global natural coordinates\n");
		output_global_parameters(&info);
		printf("Finish geometry only output\n");
		return 0;
	}
	if (isGeometryOnly)
	{
		// make connectivity
		printf("start Make_connectivity\n");
		Allocation(7, &info);
		Allocation(8, &info);
		Make_connectivity(&info);

		printf("Geometry only output mode\n\n");
		output_for_paraview_timestep(&info, isGeometryOnly, 0.0);
		printf("Finish geometry only output\n");
		return 0;
	}

	// make connectivity
	printf("start Make_connectivity\n\n");
	Allocation(7, &info);
	Allocation(8, &info);
	Make_connectivity(&info);

	printf("start output_for_paraview_timestep\n\n");
	output_for_paraview_timestep(&info, isGeometryOnly, 0.0); // start time = 0.0

	// first Inverse mapping
	printf("start first inverse mapping\n\n");
	First_inverse_mapping(&info);

	// second Inverse mapping
	printf("start second inverse mapping\n\n");
	Second_inverse_mapping(&info);

	// output as new imput file
	printf("start output new input file\n\n");
	Output_new_input_file(&info);

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\tPostprocess time:\t%.3f[s]\n\n", time);

	end[0] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[0] - start[0]).count()) / 1000.0;
	printf("\tAll analysis time:\t%.3f[s]\n", time);

	return 0;
}