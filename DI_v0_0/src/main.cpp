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
	Total_mesh = argc - 2;

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
	else if (Total_mesh >= 2)
		printf("SS-IGA carried out.(%d local meshes)\n\n", Total_mesh - 1);

	printf("start Get Input Data\n\n");

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 1st time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_1(i, argv[i + 1], &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 2nd time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_2(i, argv[i + 1], &info);

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

	if (info.c.FRACTURE_MODE == 0 && info.c.CALCLATE_DISPLACEMENT == 0)
	{
		printf("error: No valid analysis mode selected.\n");
		exit(0);
	}

	// MAX_K_WHOLE_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

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

	// make K matrix structure
	Make_D_Matrix(&info);
	Make_Index_Dof(&info);
	Make_K_Whole_Ptr_Col(&info, 0);

	// memory allocation
	Allocation(alloc_count++, &info);

	// make K matrix structure
	printf("start Make_K_Whole_Ptr_Col\n");
	Make_K_Whole_Ptr_Col(&info, 1);

	// memory allocation
	Allocation(alloc_count++, &info);

	// make connectivity
	printf("start Make_connectivity\n");
	Allocation(alloc_count++, &info);
	Allocation(alloc_count++, &info);
	Make_connectivity(&info);

	printf("start output_for_paraview_timestep\n");
	output_for_paraview_timestep(&info, isGeometryOnly, 0.0); // start time = 0.0

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\tPostprocess time:\t%.3f[s]\n\n", time);

	end[0] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[0] - start[0]).count()) / 1000.0;
	printf("\tAll analysis time:\t%.3f[s]\n", time);

	return 0;
}
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
	Total_mesh = argc - 2;

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
	else if (Total_mesh >= 2)
		printf("SS-IGA carried out.(%d local meshes)\n\n", Total_mesh - 1);

	printf("start Get Input Data\n\n");

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 1st time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_1(i, argv[i + 1], &info);

	// memory allocation
	Allocation(alloc_count++, &info);

	// Read file 2nd time
	for (int i = 0; i < Total_mesh; i++)
		Get_Input_2(i, argv[i + 1], &info);

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

	if (info.c.FRACTURE_MODE == 0 && info.c.CALCLATE_DISPLACEMENT == 0)
	{
		printf("error: No valid analysis mode selected.\n");
		exit(0);
	}

	// MAX_K_WHOLE_SIZE
	Global_var(glo_var_count++, &info);

	// memory allocation
	Allocation(alloc_count++, &info);

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

	// make K matrix structure
	Make_D_Matrix(&info);
	Make_Index_Dof(&info);
	Make_K_Whole_Ptr_Col(&info, 0);

	// memory allocation
	Allocation(alloc_count++, &info);

	// make K matrix structure
	printf("start Make_K_Whole_Ptr_Col\n");
	Make_K_Whole_Ptr_Col(&info, 1);

	// memory allocation
	Allocation(alloc_count++, &info);

	// make connectivity
	printf("start Make_connectivity\n");
	Allocation(alloc_count++, &info);
	Allocation(alloc_count++, &info);
	Make_connectivity(&info);

	printf("start output_for_paraview_timestep\n");
	output_for_paraview_timestep(&info, isGeometryOnly, 0.0); // start time = 0.0

	end[1] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[1] - start[1]).count()) / 1000.0;
	printf("\tPostprocess time:\t%.3f[s]\n\n", time);

	end[0] = chrono::system_clock::now();
	time = (double)(chrono::duration_cast<chrono::milliseconds>(end[0] - start[0]).count()) / 1000.0;
	printf("\tAll analysis time:\t%.3f[s]\n", time);

	return 0;
}