// header
#include "_header.hpp"
#include "_sub.hpp"

using namespace std;

// option file
void Get_Option(const char *s, vector<string> &opt_files)
{
	ifstream file(s);

	string line;
	while (getline(file, line))
	{
		if (line.empty() || line[0] == '#')
			continue;

		size_t commentPos = line.find('#');
		if (commentPos != string::npos)
			line = line.substr(0, commentPos);

		istringstream iss(line);
		string filename;
		if (iss >> filename)
			opt_files.emplace_back(filename);
	}
}


void Set_Option_Data(information *info)
{
	int temp_i;

	// option file 0: constant.ini is already read
	// option file 1: crack_pair.txt
	const std::string &path = info->opt_files[1];
	std::ifstream ifs(path);
	if (!ifs.is_open())
	{
		fprintf(stderr, "Error: cannot open option file '%s'\n", path.c_str());
		exit(1);
	}

	// read all integers from file, ignore comments that start with '#'
	std::vector<int> tokens;
	std::string line;
	while (std::getline(ifs, line))
	{
		// strip comment
		size_t commentPos = line.find('#');
		if (commentPos != std::string::npos)
			line = line.substr(0, commentPos);

		std::istringstream iss(line);
		int v;
		while (iss >> v)
			tokens.push_back(v);
	}

	size_t idx = 0;
	if (tokens.size() < static_cast<size_t>(Total_mesh))
	{
		fprintf(stderr, "Error: not enough data in '%s' (need %d mesh counts)\n", path.c_str(), Total_mesh);
		exit(1);
	}

	info->crack_pair_n_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));
	if (info->crack_pair_n_to_mesh == NULL)
	{
		fprintf(stderr, "Cannot allocate memory for crack_pair_n_to_mesh\n");
		exit(1);
	}

	int total_crack_pair_n = 0;
	for (int i = 0; i < Total_mesh; i++)
	{
		temp_i = tokens[idx++];
		info->crack_pair_n_to_mesh[i] = total_crack_pair_n;
		// guard against negative counts
		if (temp_i < 0)
		{
			fprintf(stderr, "Error: negative crack pair count in mesh %d\n", i);
			exit(1);
		}
		total_crack_pair_n += temp_i;
	}
	info->crack_pair_n_to_mesh[Total_mesh] = total_crack_pair_n;

	// need total_crack_pair_n * 2 integers remaining
	if (tokens.size() - idx < static_cast<size_t>(total_crack_pair_n) * 2)
	{
		fprintf(stderr, "Error: not enough crack pair entries in '%s' (expected %d pairs)\n", path.c_str(), total_crack_pair_n);
		exit(1);
	}

	info->crack_pair = (int *)malloc(sizeof(int) * total_crack_pair_n * 2);
	if (info->crack_pair == NULL && total_crack_pair_n > 0)
	{
		fprintf(stderr, "Cannot allocate memory for crack_pair\n");
		exit(1);
	}

	for (int i = 0; i < total_crack_pair_n; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			temp_i = tokens[idx++];
			info->crack_pair[i * 2 + j] = temp_i;
		}
	}
}


// memory allocation
void Allocation(const int num, information *info)
{
	if (num == 0)
	{
		// general
		info->c.GEOMETRY_ONLY_OUTPUT          = any_cast<int>(info->c.data_vec[0]);
		info->c.OUTPUT_GLOBAL_PARAMETERS      = any_cast<int>(info->c.data_vec[1]);
		info->c.CALC_ON_ELE_VERTEX            = any_cast<int>(info->c.data_vec[2]);

		// gaussian quadrature
		info->c.USE_EXTENDED_QUADRATURE       = any_cast<int>(info->c.data_vec[3]);
		info->c.NUM_GAUSS_POINTS              = any_cast<int>(info->c.data_vec[4]);
		info->c.NUM_GAUSS_POINTS_EXTENDED     = any_cast<int>(info->c.data_vec[5]);

		info->Total_Knot_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		info->Total_Patch_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));			  // 各メッシュ上のパッチ数
		info->Total_Patch_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		  // メッシュまでのパッチ数（メッシュ内のパッチ数は含まない）
		info->Total_Control_Point_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));	  // 各メッシュ上のコントロールポイント数
		info->Total_Control_Point_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int)); // メッシュまでのコントロールポイント数(メッシュ内のコントロールポイント数は含まない)
		info->Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		  // 各メッシュ上の要素数
		info->Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		  // メッシュまでの要素数(メッシュ内の要素数は含まない)
		info->real_Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));	  // real_Total_Element_on_mesh[MAX_N_MESH]
		info->real_Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));  // real_Total_Element_to_mesh[MAX_N_MESH + 1]
		info->Total_Load_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		info->Total_Constraint_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		info->Total_DistributeForce_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
		if (info->Total_Knot_to_mesh == NULL || info->Total_Patch_on_mesh == NULL || info->Total_Patch_to_mesh == NULL || info->Total_Control_Point_on_mesh == NULL || info->Total_Control_Point_to_mesh == NULL || info->Total_Element_on_mesh == NULL || info->Total_Element_to_mesh == NULL || info->real_Total_Element_on_mesh == NULL || info->real_Total_Element_to_mesh == NULL || info->Total_Load_to_mesh == NULL || info->Total_Constraint_to_mesh == NULL || info->Total_DistributeForce_to_mesh == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}

		int max_pow_ng = 0;
		if (info->c.USE_EXTENDED_QUADRATURE == 0)
			max_pow_ng = pow_int(info->c.NUM_GAUSS_POINTS, MAX_DIMENSION);
		else if (info->c.USE_EXTENDED_QUADRATURE == 1)
			max_pow_ng = pow_int(info->c.NUM_GAUSS_POINTS_EXTENDED, MAX_DIMENSION);
		else
		{
			printf("info->c.USE_EXTENDED_QUADRATURE is not defined\n");
			exit(1);
		}

		info->gauss_w = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point = (double *)malloc(sizeof(double) * max_pow_ng * MAX_DIMENSION);
		info->gauss_w_ex = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point_ex = (double *)malloc(sizeof(double) * max_pow_ng * MAX_DIMENSION);
		info->gauss_w_1D = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point_1D = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_w_1D_ex = (double *)malloc(sizeof(double) * max_pow_ng);
		info->gauss_point_1D_ex = (double *)malloc(sizeof(double) * max_pow_ng);
		if (info->gauss_w == NULL || info->gauss_point == NULL || info->gauss_w_ex == NULL || info->gauss_point_ex == NULL || info->gauss_w_1D == NULL || info->gauss_point_1D == NULL || info->gauss_w_1D_ex == NULL || info->gauss_point_1D_ex == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 1)
	{
		info->Order = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));								  // Order[MAX_N_PATCH][info->DIMENSION]
		info->No_knot = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));							  // No_knot[MAX_N_PATCH][info->DIMENSION]
		info->Total_Control_Point_to_patch = (int *)calloc((info->Total_Patch_to_mesh[Total_mesh] + 1), sizeof(int));					  // Total_Control_Point_to_patch[MAX_N_PATCH]
		info->Total_Knot_to_patch_dim = (int *)calloc((info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION + 1), sizeof(int));		  // Total_Knot_to_patch_dim[MAX_N_PATCH][info->DIMENSION]
		info->Position_Knots = (double *)malloc(sizeof(double) * info->Total_Knot_to_mesh[Total_mesh]);									  // Position_Knots[MAX_N_PATCH][info->DIMENSION][MAX_N_KNOT];
		info->No_Control_point = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));					  // Order[MAX_N_PATCH][info->DIMENSION]
		info->No_Control_point_in_patch = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh]));							  // No_Control_point_in_patch[MAX_N_PATCH]
		info->Patch_Control_point = (int *)malloc(sizeof(int) * MAX_CP);																  // Patch_Control_point[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]
		info->No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh]));						  // No_Control_point_ON_ELEMENT[MAX_N_PATCH]
		info->Node_Coordinate = (double *)malloc(sizeof(double) * (info->Total_Control_Point_to_mesh[Total_mesh] * (info->DIMENSION + 1))); // Node_Coordinate[MAX_N_NODE][info->DIMENSION + 1];
		info->Control_Coord_x = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Coord_y = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Coord_z = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Weight = (double *)malloc(sizeof(double) * MAX_CP);																  // Control_Weight[MAX_N_NODE];
		info->Constraint_Node_Dir = (int *)malloc(sizeof(int) * (info->Total_Constraint_to_mesh[Total_mesh] * 2));						  // Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
		info->Value_of_Constraint = (double *)calloc(info->Total_Constraint_to_mesh[Total_mesh], sizeof(double));						  // Value_of_Constraint[MAX_N_CONSTRAINT];
		info->Load_Node_Dir = (int *)malloc(sizeof(int) * (info->Total_Load_to_mesh[Total_mesh] * 2));									  // Load_Node_Dir[MAX_N_LOAD][2];
		info->Value_of_Load = (double *)calloc(info->Total_Load_to_mesh[Total_mesh], sizeof(double));									  // Value_of_Load[MAX_N_LOAD];
		info->iPatch_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // iPatch_array[MAX_N_DISTRIBUTE_FORCE]
		info->iCoord_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // iCoord_array[MAX_N_DISTRIBUTE_FORCE]
		info->jCoord_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);								  // jCoord_array[MAX_N_DISTRIBUTE_FORCE]
		info->type_load_array = (int *)malloc(sizeof(int) * info->Total_DistributeForce_to_mesh[Total_mesh]);							  // type_load_array[MAX_N_DISTRIBUTE_FORCE]
		info->val_Coord_array = (double *)calloc(info->Total_DistributeForce_to_mesh[Total_mesh], sizeof(double));						  // val_Coord_array[MAX_N_DISTRIBUTE_FORCE]
		info->Range_Coord_array = (double *)calloc((info->Total_DistributeForce_to_mesh[Total_mesh] * 2 * 2), sizeof(double));			  // Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2 * 2] (info->DIMENSION == 3 のとき 2 -> 4)
		info->Coeff_Dist_Load_array = (double *)calloc((info->Total_DistributeForce_to_mesh[Total_mesh] * 3 * 2), sizeof(double));		  // Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3 * 2] (info->DIMENSION == 3 のとき 3 -> 6)
		if (info->Order == NULL || info->No_knot == NULL || info->Total_Control_Point_to_patch == NULL || info->Total_Knot_to_patch_dim == NULL || info->Position_Knots == NULL || info->No_Control_point == NULL || info->No_Control_point_in_patch == NULL || info->Patch_Control_point == NULL || info->No_Control_point_ON_ELEMENT == NULL || info->Node_Coordinate == NULL || info->Control_Coord_x == NULL || info->Control_Coord_y == NULL || info->Control_Weight == NULL || info->Constraint_Node_Dir == NULL || info->Value_of_Constraint == NULL || info->Load_Node_Dir == NULL || info->Value_of_Load == NULL || info->iPatch_array == NULL || info->iCoord_array == NULL || info->jCoord_array == NULL || info->type_load_array == NULL || info->val_Coord_array == NULL || info->Range_Coord_array == NULL || info->Coeff_Dist_Load_array == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
		info->Geo_Order = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->DIMENSION));								  // Order[MAX_N_PATCH][info->DIMENSION]
		info->Geo_No_knot = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->DIMENSION));							  // No_knot[MAX_N_PATCH][info->DIMENSION]
		info->Geo_Total_Control_Point_to_patch = (int *)calloc((info->Geo_Total_patch_on_mesh + 1), sizeof(int));					  // Total_Control_Point_to_patch[MAX_N_PATCH]
		info->Geo_Total_Knot_to_patch_dim = (int *)calloc((info->Geo_Total_patch_on_mesh * info->DIMENSION + 1), sizeof(int));
		info->Geo_Position_Knots = (double *)malloc(sizeof(double) * info->Geo_Total_Knot_on_mesh);									  // Position_Knots[MAX_N_PATCH][info->DIMENSION][MAX_N_KNOT];
		info->Geo_No_Control_point = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->DIMENSION));
		info->Geo_No_Control_point_in_patch = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh));
		info->Geo_Patch_Control_point = (int *)malloc(sizeof(int) * MAX_CP);
		info->Geo_No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh));
		info->Geo_Node_Coordinate = (double *)malloc(sizeof(double) * (info->Geo_Total_Control_Point_on_mesh * (info->DIMENSION + 1))); 	// Geo_Node_Coordinate[MAX_N_NODE][info->DIMENSION + 1];
		info->Control_Coord_x_geo = (double *)malloc(sizeof(double) * MAX_CP);																  			// Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Coord_y_geo = (double *)malloc(sizeof(double) * MAX_CP);			
		info->Control_Coord_z_geo = (double *)malloc(sizeof(double) * MAX_CP);																  			// Control_Coord[info->DIMENSION][MAX_N_NODE];
		info->Control_Weight_geo = (double *)malloc(sizeof(double) * MAX_CP);																  			// Control_Weight[MAX_N_NODE];
		if (info->Geo_Order == NULL || info->Geo_No_knot == NULL || info->Geo_Total_Control_Point_to_patch == NULL || info->Geo_Total_Knot_to_patch_dim == NULL || info->Geo_Position_Knots == NULL || info->Geo_No_Control_point == NULL || info->Geo_No_Control_point_in_patch == NULL || info->Geo_Patch_Control_point == NULL || info->Geo_No_Control_point_ON_ELEMENT == NULL || info->Geo_Node_Coordinate == NULL || info->Control_Coord_x_geo == NULL || info->Control_Coord_y_geo == NULL || info->Control_Weight_geo == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 2)
	{
		info->INC = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION));		 // INC[MAX_N_PATCH][MAX_N_NODE][info->DIMENSION]
		info->Controlpoint_of_Element = (int *)malloc(sizeof(int) * (info->Total_Element_to_mesh[Total_mesh] * MAX_NO_CP_ON_ELEMENT));							 // Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CP_ON_ELEMENT]
		info->Element_patch = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																 // Element_patch[MAX_N_ELEMENT]
		info->Element_mesh = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																 // Element_mesh[MAX_N_ELEMENT] 要素がどのメッシュ内にあるかを示す配列
		info->line_No_real_element = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));										 // line_No_real_element[MAX_N_PATCH][info->DIMENSION] ゼロエレメントではない要素列の数
		info->line_No_Total_element = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->DIMENSION));										 // line_No_Total_element[MAX_N_PATCH][info->DIMENSION] ゼロエレメントを含むすべての要素列の数
		info->difference = (double *)calloc(info->Total_Knot_to_mesh[Total_mesh], sizeof(double));																 // difference[MAX_N_PATCH][info->DIMENSION][MAX_N_KNOT] 隣り合うノットベクトルの差
		info->Total_element_all_ID = (int *)calloc(info->Total_Element_to_mesh[Total_mesh], sizeof(int));														 // Total_element_all_ID[MAX_N_ELEMENT] ゼロエレメントではない要素 = 1, ゼロエレメント = 0
		info->ENC = (int *)malloc(sizeof(int) * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION));														 // ENC[MAX_N_ELEMENT][info->DIMENSION] ENC[パッチ][全ての要素][0, 1] = x, y方向の何番目の要素か
		info->real_element_line = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION)); // real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][info->DIMENSION] ゼロエレメントではない要素列
		info->real_element = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh]);																 // real_element[MAX_N_ELEMENT] ゼロエレメントではない要素の番号
		info->real_El_No_on_mesh = (int *)malloc(sizeof(int) * (info->Total_Patch_to_mesh[Total_mesh] * info->Total_Element_to_mesh[Total_mesh]));				 // real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT]
		info->Equivalent_Nodal_Force = (double *)calloc(MAX_CP * info->DIMENSION, sizeof(double));																 // Equivalent_Nodal_Force[MAX_N_NODE][info->DIMENSION] Equivalent nodal forces arising from the distributed load
		if (info->INC == NULL || info->Controlpoint_of_Element == NULL || info->Element_patch == NULL || info->Element_mesh == NULL || info->line_No_real_element == NULL || info->line_No_Total_element == NULL || info->difference == NULL || info->Total_element_all_ID == NULL || info->ENC == NULL || info->real_element_line == NULL || info->real_element == NULL || info->real_El_No_on_mesh == NULL || info->Equivalent_Nodal_Force == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
		if (Total_mesh >= 2)
		{
			info->Geo_INC = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->Geo_Total_Control_Point_to_patch[info->Geo_Total_patch_on_mesh] * info->DIMENSION));
			info->Geo_Controlpoint_of_Element = (int *)malloc(sizeof(int) * (info->Geo_Total_Element_on_mesh * MAX_NO_CP_ON_ELEMENT));
			info->Geo_Element_patch = (int *)malloc(sizeof(int) * info->Geo_Total_Element_on_mesh);
			info->Geo_Element_mesh = (int *)malloc(sizeof(int) * info->Geo_Total_Element_on_mesh);
			info->Geo_line_No_real_element = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->Geo_Total_Element_on_mesh * info->DIMENSION));
			info->Geo_line_No_Total_element = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->Geo_Total_Element_on_mesh * info->DIMENSION));
			info->Geo_difference = (double *)calloc(info->Geo_Total_Knot_on_mesh, sizeof(double));
			info->Geo_Total_element_all_ID = (int *)calloc(info->Geo_Total_Element_on_mesh, sizeof(int));
			info->Geo_ENC = (int *)malloc(sizeof(int) * (info->Geo_Total_Element_on_mesh * info->DIMENSION));
			info->Geo_real_element_line = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->Geo_Total_Element_on_mesh * info->DIMENSION));
			info->Geo_real_element = (int *)malloc(sizeof(int) * info->Geo_Total_Element_on_mesh);
			info->Geo_real_El_No_on_mesh = (int *)malloc(sizeof(int) * (info->Geo_Total_patch_on_mesh * info->Geo_Total_Element_on_mesh));
			info->Geo_Equivalent_Nodal_Force = (double *)calloc(MAX_CP * info->DIMENSION, sizeof(double));
			if (info->Geo_INC == NULL || info->Geo_Controlpoint_of_Element == NULL || info->Geo_Element_patch == NULL || info->Geo_Element_mesh == NULL || info->Geo_line_No_real_element == NULL || info->Geo_line_No_Total_element == NULL || info->Geo_difference == NULL || info->Geo_Total_element_all_ID == NULL || info->Geo_ENC == NULL || info->Geo_real_element_line == NULL || info->Geo_real_element == NULL || info->Geo_real_El_No_on_mesh == NULL || info->Geo_Equivalent_Nodal_Force == NULL)
			{
				printf("Cannot allocate memory\n");
				exit(1);
			}
		}
	}
	else if (num == 3)
	{
		info->gp.resize(info->Total_Element_to_mesh[Total_mesh]);
		info->eoi.resize(info->Total_Element_to_mesh[Total_mesh]);
	}
	else if (num == 4)
	{
		info->D = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * D_MATRIX_SIZE);
		info->Index_Dof = (int *)calloc(MAX_K_WHOLE_SIZE, sizeof(int));													   // Index_Dof[MAX_K_WHOLE_SIZE];
		info->K_Whole_Ptr = (long long *)calloc(MAX_K_WHOLE_SIZE + 1, sizeof(long long));										   // K_Whole_Ptr[MAX_K_WHOLE_SIZE + 1]
		if (info->D == NULL || info->Index_Dof == NULL || info->K_Whole_Ptr == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 5)
	{
		info->K_Whole_Col = (int *)malloc(sizeof(int) * info->K_Whole_Ptr[K_Whole_Size]);	  // K_Whole_Col[MAX_NON_ZERO]
		info->K_Whole_Val = (double *)calloc(info->K_Whole_Ptr[K_Whole_Size], sizeof(double)); // K_Whole_Val[MAX_NON_ZERO]
		if (info->K_Whole_Col == NULL || info->K_Whole_Val == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 6)
	{
		info->sol_vec = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(double)); // sol_vec[MAX_K_WHOLE_SIZE]
		info->rhs_vec = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(double)); // rhs_vec[MAX_K_WHOLE_SIZE]
		if (info->sol_vec == NULL || info->rhs_vec == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 7)
	{
		info->Displacement = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE);
		if (info->Displacement == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if (num == 8)
	{
		int total_array = 0, check_n = 0, FE_n = 0, opp_n = 0, point_on_element = pow_int(3, info->DIMENSION);
		int total_edge = 0, edge_on_ele = 0;
		if (info->DIMENSION == 2)
		{
			for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
			{
				int x = info->line_No_Total_element[i * info->DIMENSION] * 2 + 1, y = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1;
				total_array += 2 * (x + y);
			}
			check_n = 16;
			FE_n = 32;
			opp_n = 4;
			edge_on_ele = 2;
		}
		else if (info->DIMENSION == 3)
		{
			for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
			{
				int x = info->line_No_Total_element[i * info->DIMENSION] * 2 + 1, y = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1, z = info->line_No_Total_element[i * info->DIMENSION + 2] * 2 + 1;
				total_array += 2 * (x * y + y * z + x * z);
			}
			check_n = 24;
			FE_n = 36;
			opp_n = 6;
			edge_on_ele = 4;
		}
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				total_edge += info->line_No_Total_element[i * info->DIMENSION + j] * edge_on_ele * 3;
			}
		}
		info->Connectivity = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->Connectivity_ele = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->Connectivity_point = (int *)malloc(sizeof(int) * info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->Connectivity_coord = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * info->DIMENSION, sizeof(double));
		info->Patch_check = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * check_n);
		info->Patch_array = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * total_array);
		info->Face_Edge_info = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * FE_n);
		info->Opponent_patch_num = (int *)malloc(sizeof(int) * info->Total_Patch_to_mesh[Total_mesh] * opp_n);
		info->Edge_coord = (double *)calloc(total_edge * info->DIMENSION, sizeof(double));
		info->disp_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * info->DIMENSION, sizeof(double));
		info->strain_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * N_STRAIN, sizeof(double));
		info->stress_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element * N_STRESS, sizeof(double));
		info->vm_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element, sizeof(double));
		info->hs_at_connectivity = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * point_on_element, sizeof(double));
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh] * point_on_element; i++)
		{
			info->Connectivity[i] = -1;
		}
		for (int i = 0; i < info->Total_Patch_to_mesh[Total_mesh] * FE_n; i++)
		{
			info->Face_Edge_info[i] = -1;
		}
		if (info->Connectivity == NULL || info->Connectivity_ele == NULL || info->Connectivity_point == NULL || info->Connectivity_coord == NULL || info->Patch_check == NULL || info->Patch_array == NULL || info->Face_Edge_info == NULL || info->Opponent_patch_num == NULL || info->disp_at_connectivity == NULL || info->strain_at_connectivity == NULL || info->stress_at_connectivity == NULL || info->vm_at_connectivity == NULL || info->hs_at_connectivity == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}

		info->Connectivity_all_ele.resize(info->Total_Element_to_mesh[Total_mesh] * point_on_element);
		info->ecn.resize(info->Total_Element_to_mesh[Total_mesh]);
	}
	else if(num == 9)
	{
		int n_gp = pow_int(info->c.NUM_GAUSS_POINTS, info->DIMENSION);
		int n_gp_ex = pow_int(info->c.NUM_GAUSS_POINTS_EXTENDED, info->DIMENSION);
		int max_gp_per_element = (n_gp_ex > n_gp) ? n_gp_ex : n_gp;
		info->Displacement_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * max_gp_per_element * info->DIMENSION, sizeof(double));
		info->PhysicalCoordinate_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * max_gp_per_element * info->DIMENSION, sizeof(double));
		info->Strain_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * max_gp_per_element * N_STRAIN, sizeof(double));
		info->Stress_at_GP = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * max_gp_per_element * N_STRESS, sizeof(double));
		info->Displacement_at_ele_vertex = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * pow_int(3, info->DIMENSION) * info->DIMENSION, sizeof(double));
		info->Strain_at_ele_vertex = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * pow_int(3, info->DIMENSION) * N_STRAIN, sizeof(double));
		info->Stress_at_ele_vertex = (double *)calloc(info->Total_Element_to_mesh[Total_mesh] * pow_int(3, info->DIMENSION) * N_STRESS, sizeof(double));
		info->PhysicalCoordinate_at_ele_vertex = (double *)calloc(info->Geo_Total_Element_on_mesh * pow_int(3, info->DIMENSION) * info->DIMENSION, sizeof(double));
		info->ReactionForce = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(double)); // ReactionForce[MAX_K_WHOLE_SIZE]
		if (info->Strain_at_GP == NULL || info->Stress_at_GP == NULL || info->ReactionForce == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else if(num == 10)
	{
		// memory allocation
		int Element_vertex = info->Geo_Total_Element_on_mesh;
		Element_vertex *= pow_int(3, info->DIMENSION);
		info->Target_Physical_Coord = (double *)malloc(sizeof(double) * Element_vertex * info->DIMENSION); // Target_Physical_Coord[Element_vertex][info->DIMENSION]
		info->Target_Para_Coord = (double *)malloc(sizeof(double) * Element_vertex * info->DIMENSION); // Target_Para_Coord[Element_vertex][info->DIMENSION]
		info->New_Node_Coordinate = (double *)malloc(sizeof(double) * (info->Geo_Total_Control_Point_on_mesh * (info->DIMENSION + 1)));
		if (info->Target_Physical_Coord == NULL || info->Target_Para_Coord == NULL || info->New_Node_Coordinate == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(1);
		}
	}
	else
	{
		printf("Invalid num in Global_var function\n");
		exit(1);
	}
}

// set global variable
void Global_var(const int num, information *info)
{
	if (num == 0)
	{
		MAX_ORDER += 1; // 次数の最大値 + 1
		MAX_NO_CP_ON_ELEMENT = pow_int(MAX_ORDER, info->DIMENSION);
		MAX_KIEL_SIZE = MAX_NO_CP_ON_ELEMENT * info->DIMENSION;
	}
	else if (num == 1)
	{
		if (info->DIMENSION == 2)
		{
			D_MATRIX_SIZE = 3;

			// 0: xx, 1: yy, 2: xy, 3: zz
			N_STRAIN = 4;
			N_STRESS = 4;
		}
		else if (info->DIMENSION == 3)
		{
			D_MATRIX_SIZE = 6;

			// 0: xx, 1: yy, 2: zz, 3: xy, 4: yz, 5: xz
			N_STRAIN = 6;
			N_STRESS = 6;
		}
	}
	else if (num == 2)
	{
		MAX_K_WHOLE_SIZE = info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION;
	}
}


// Read file 1st time
void Get_Input_1(int tm, const char *filename, information *info)
{
	char s[256];
	int temp_i;

	int i, j;
	
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Error: Cannot open input file %s\n", filename);
		exit(1);
	}

	printf("\n----- Reading Input File %d -----\n", tm + 1);

	// DIMENSION
	fscanf(fp, "%d", &temp_i);
	info->DIMENSION = temp_i;
	fgets(s, 256, fp);
	printf("DIMENSION: %d\n", info->DIMENSION);
	if (info->DIMENSION != 2 && info->DIMENSION != 3)
	{
		printf("Error, wrong DIMENSION in input data\n");
		exit(1);
	}

	// 材料定数
	fscanf(fp, "%lf %lf", &E, &nu);
	fgets(s, 256, fp);
	printf("E: %le, nu: %le\n", E, nu);

	// パッチ数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;
	int *CP = (int *)malloc(sizeof(int) * temp_i * info->DIMENSION);
	if (CP == NULL) {
		printf("Error: Memory allocation failed for CP array\n");
		fclose(fp);
		exit(1);
	}
	printf("No_Patch: %d\n", temp_i);
	info->Total_Patch_on_mesh[tm] = temp_i;
	info->Total_Patch_to_mesh[tm + 1] = info->Total_Patch_to_mesh[tm] + temp_i;
	printf("Total_Patch_to_mesh[%d] = %d\n", tm, info->Total_Patch_to_mesh[tm]);

	// コントロールポイント数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	printf("Total_Control_Point: %d\n", temp_i);
	info->Total_Control_Point_on_mesh[tm] = temp_i;
	info->Total_Control_Point_to_mesh[tm + 1] = info->Total_Control_Point_to_mesh[tm] + temp_i;
	printf("Total_Control_Point_to_mesh[%d] = %d\n", tm, info->Total_Control_Point_to_mesh[tm]);

	// 各方向の次数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	info->Total_Knot_to_mesh[tm + 1] = info->Total_Knot_to_mesh[tm];
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Total_Knot_to_mesh[tm + 1] += temp_i;
		}
	}

	// 各パッチ各方向のコントロールポイント数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			CP[i * info->DIMENSION + j] = temp_i;
		}

		int temp_MAX_CP = 1;
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_MAX_CP *= CP[i * info->DIMENSION + j];
		}
		MAX_CP += temp_MAX_CP;
	}
	fgets(s, 256, fp);

	// パッチコネクティビティ(スキップ)
	for (i = 0; i < info->Total_Patch_on_mesh[tm]; i++)
	{
		int temp_MAX_CP = 1;
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_MAX_CP *= CP[i * info->DIMENSION + j];
		}
		for (j = 0; j < temp_MAX_CP; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// 境界条件
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	info->Total_Constraint_to_mesh[tm + 1] = info->Total_Constraint_to_mesh[tm] + Total_Constraint;
	info->Total_Load_to_mesh[tm + 1] = info->Total_Load_to_mesh[tm] + Total_Load;
	info->Total_DistributeForce_to_mesh[tm + 1] = info->Total_DistributeForce_to_mesh[tm] + Total_DistributeForce;
	
	printf("Total_Constraint = %d\n", Total_Constraint);
	printf("Total_Load = %d\n", Total_Load);
	printf("Total_DistributedForce = %d\n", Total_DistributeForce);
	fgets(s, 256, fp);
	
	if (tm == 1)
	{
		// ノットベクトルの読み込み	(スキップ)
		for (i = 0; i < No_Patch; i++)
		{
		    for (j = 0; j < info->DIMENSION; j++)
		    {			
		        int c;
		        while ((c = fgetc(fp)) != '\n' && c != EOF) {
		        }
		    }
		}
		int c;
		while ((c = fgetc(fp)) != '\n' && c != EOF) {
		}

		// 境界条件の詳細(スキップ)
		for (i = 0; i < Total_Constraint; i++)
		{
			int c;
		    while ((c = fgetc(fp)) != '\n' && c != EOF) {
		    }
		}
		fgets(s, 256, fp);

		// 境界条件の詳細(スキップ)
		for (i = 0; i < Total_DistributeForce; i++)
		{
			int c;
		    while ((c = fgetc(fp)) != '\n' && c != EOF) {
		    }
		}
		fgets(s, 256, fp);

		// ローカルパッチが属するグローバルパッチ番号
		fscanf(fp, "%d", &temp_i);
		info->Global_local_patch = temp_i;
		printf("Global patch number, the local patch belongs to = %d\n", info->Global_local_patch);
		fgets(s, 256, fp);

		// ジオメトリを定義するパッチ数
		fscanf(fp, "%d", &temp_i);
		info->Geo_Total_patch_on_mesh = temp_i;
		No_Patch = temp_i;
		printf("Geometry patch on mesh = %d\n", info->Geo_Total_patch_on_mesh);

		// ジオメトリを定義する制御点の数
		fscanf(fp, "%d", &temp_i);
		info->Geo_Total_Control_Point_on_mesh = temp_i;
		printf("Geometry Control point on mesh = %d\n", info->Geo_Total_Control_Point_on_mesh);

		// 各方向の次数(スキップ)
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				fscanf(fp, "%d", &temp_i);
			}
		}
		fgets(s, 256, fp);

		// ノット数 (スキップ)
		info->Geo_Total_Knot_on_mesh = 0;
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				fscanf(fp, "%d", &temp_i);
				info->Geo_Total_Knot_on_mesh += temp_i;
			}
		}

		// 各パッチ各方向のコントロールポイント数(スキップ)
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				fscanf(fp, "%d", &temp_i);
				CP[i * info->DIMENSION + j] = temp_i;
			}

			int temp_MAX_CP = 1;
			for (j = 0; j < info->DIMENSION; j++)
			{
				temp_MAX_CP *= CP[i * info->DIMENSION + j];
			}
			MAX_CP += temp_MAX_CP;
		}
		fgets(s, 256, fp);


		// パッチコネクティビティ(スキップ)
		for (i = 0; i < info->Total_Patch_on_mesh[tm]; i++)
		{
			int temp_MAX_CP = 1;
			for (j = 0; j < info->DIMENSION; j++)
			{
				temp_MAX_CP *= CP[i * info->DIMENSION + j];
			}
			for (j = 0; j < temp_MAX_CP; j++)
			{
				fscanf(fp, "%d", &temp_i);
			}
		}
		fgets(s, 256, fp);

		// 境界条件 (スキップ)
		int Total_Constraint, Total_Load, Total_DistributeForce;
		fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
		fgets(s, 256, fp);

		printf("\n");

		}
		
	fclose(fp);
	free(CP);
}


// Read file 2nd time
void Get_Input_2(int tm, const char *filename, information *info)
{
	char s[256];
	int temp_i;
	double temp_d;

	int i, j, k;

	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Error: Cannot open input file %s\n", filename);
		exit(1);
	}

	// DIMENSION(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);

	// 材料定数(スキップ)
	fscanf(fp, "%lf %lf", &temp_d, &temp_d);
	fgets(s, 256, fp);

	// パッチ数(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;

	// コントロールポイント数(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int Total_Control_Point = temp_i;

	// 各方向の次数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			if (MAX_ORDER < temp_i)
			{
				MAX_ORDER = temp_i;
			}
			printf("Order[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			if (MAX_KNOT < temp_i)
			{
				MAX_KNOT = temp_i;
			}
			info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j + 1] = info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + temp_i;
			printf("No_knot[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// 各パッチ各方向のコントロールポイント数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			printf("No_Control_point[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	for (i = 0; i < No_Patch; i++)
	{
		info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]] = 1;
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]] *= info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j];
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm] + 1] = info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]];
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			if (info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] != info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in mesh_No.%d in patch_No.%d direction:%d\n", tm, i, j);
			}
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_in_patch[%d] = %d\t", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]);
	}
	printf("\n");

	// パッチコネクティビティ
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Patch_Control_point[info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + j] = temp_i;
			if (tm > 0)
			{
				info->Patch_Control_point[info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + j] += info->Total_Control_Point_to_mesh[tm];
			}
		}
	}
	fgets(s, 256, fp);

	// 境界条件(スキップ)
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	fgets(s, 256, fp);

	// ノットベクトルの読み込み
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]; k++)
			{
				fscanf(fp, "%lf", &temp_d);
				info->Position_Knots[info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + k] = temp_d;
			}
		}
	}
	fgets(s, 256, fp);

	int Total_Element = 0;
	for (i = 0; i < No_Patch; i++)
	{
		if (info->DIMENSION == 2)
		{
			Total_Element += (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1]);
			info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]] = (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] + 1);
		}
		else if (info->DIMENSION == 3)
		{
			Total_Element += (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2]);
			info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]] = (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2] + 1);
		}
	}
	printf("Total_Element = %d\n", Total_Element);
	info->Total_Element_on_mesh[tm] = Total_Element;
	info->Total_Element_to_mesh[tm + 1] = info->Total_Element_to_mesh[tm] + Total_Element;
	printf("Total_Element_on_mesh[%d] = %d\n", tm, info->Total_Element_on_mesh[tm]);

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_ON_ELEMENT[%d] = %d\n", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]]);
	}


	// 節点座標
	if (tm == 0)
	{
		for (i = 0; i < Total_Control_Point; i++)
		{
		    fscanf(fp, "%d", &temp_i);
		    for (j = 0; j < info->DIMENSION + 1; j++)
		    {
		        fscanf(fp, "%lf", &temp_d);
		        double coordinate_value = temp_d;
				info->Node_Coordinate[(temp_i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j] = coordinate_value;

			}
		}

		for (i = 0; i < No_Patch; i++)
		{
			int cp_to_patch = info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]];
			for (j = 0; j < info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]; j++)
			{
				int connectivity = info->Patch_Control_point[cp_to_patch + j] * (info->DIMENSION + 1);
				if (info->DIMENSION == 2)
				{
					info->Control_Coord_x[cp_to_patch + j] = info->Node_Coordinate[connectivity];
					info->Control_Coord_y[cp_to_patch + j] = info->Node_Coordinate[connectivity + 1];
					info->Control_Weight[cp_to_patch + j] = info->Node_Coordinate[connectivity + 2];
				}
				else if (info->DIMENSION == 3)
				{
					info->Control_Coord_x[cp_to_patch + j] = info->Node_Coordinate[connectivity];
					info->Control_Coord_y[cp_to_patch + j] = info->Node_Coordinate[connectivity + 1];
					info->Control_Coord_z[cp_to_patch + j] = info->Node_Coordinate[connectivity + 2];
					info->Control_Weight[cp_to_patch + j] = info->Node_Coordinate[connectivity + 3];
				}
			}
		}
		fgets(s, 256, fp);

		// 拘束
		for (i = 0; i < Total_Constraint; i++)
		{
			fscanf(fp, "%d %d %lf",
				   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0],
				   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 1],
				   &info->Value_of_Constraint[i + info->Total_Constraint_to_mesh[tm]]);
			info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] = info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];
		}
		if (Total_Constraint > 0)
		fgets(s, 256, fp);

		// 荷重
		for (i = 0; i < Total_Load; i++)
		{
			fscanf(fp, "%d %d %lf",
				   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
				   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
				   &info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
			info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] = info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];

			printf("Load_Node_Dir[%d]= %d Load_Node_Dir[%d]= %d Value_of_Load[%d] = %le\n",
				   (i + info->Total_Load_to_mesh[tm]) * 2 + 0, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
				   (i + info->Total_Load_to_mesh[tm]) * 2 + 1, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
				   i + info->Total_Load_to_mesh[tm], info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
		}
		if (Total_Load > 0)
		fgets(s, 256, fp);

		if (info->DIMENSION == 2)
		{
			int type_load, iPatch, iCoord;
			double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];

			for (i = 0; i < Total_DistributeForce; i++)
			{
				fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
				printf("Distibuted load number: %d\n", i);
				printf("type_load: %d iPatch: %d iCoord: %d val_Coord: %le Range_Coord: %le %le\nCoef_Dist_Load: %le %le %le\n",
						type_load, iPatch, iCoord, val_Coord, Range_Coord[0], Range_Coord[1], Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2]);

				// for S-IGA
				info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
				info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
				info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
				info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0] = Range_Coord[0];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1] = Range_Coord[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0] = Coeff_Dist_Load[0];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1] = Coeff_Dist_Load[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2] = Coeff_Dist_Load[2];
			}
		}
		else if (info->DIMENSION == 3)
		{
			int type_load, iPatch, iCoord, jCoord;
			double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];

			for (i = 0; i < Total_DistributeForce; i++)
			{
				fscanf(fp, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &jCoord, &val_Coord, &iRange_Coord[0], &iRange_Coord[1], &jRange_Coord[0], &jRange_Coord[1], &iCoeff_Dist_Load[0], &iCoeff_Dist_Load[1], &iCoeff_Dist_Load[2], &jCoeff_Dist_Load[0], &jCoeff_Dist_Load[1], &jCoeff_Dist_Load[2]);
				printf("Distibuted load number: %d\n", i);
				printf("type_load: %d iPatch: %d iCoord: %d jCoord: %d val_Coord: %le\nRange_Coord: %le %le %le %le\nCoef_Dist_Load: %le %le %le %le %le %le\n",
						type_load, iPatch, iCoord, jCoord, val_Coord,
						iRange_Coord[0], iRange_Coord[1], jRange_Coord[0], jRange_Coord[1],
						iCoeff_Dist_Load[0], iCoeff_Dist_Load[1], iCoeff_Dist_Load[2], jCoeff_Dist_Load[0], jCoeff_Dist_Load[1], jCoeff_Dist_Load[2]);

				// for S-IGA
				info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
				info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
				info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
				info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = jCoord;
				info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0] = iRange_Coord[0];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1] = iRange_Coord[1];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2] = jRange_Coord[0];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3] = jRange_Coord[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0] = iCoeff_Dist_Load[0];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1] = iCoeff_Dist_Load[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2] = iCoeff_Dist_Load[2];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3] = jCoeff_Dist_Load[0];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4] = jCoeff_Dist_Load[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5] = jCoeff_Dist_Load[2];
			}
		}
		if (Total_DistributeForce > 0)
		fgets(s, 256, fp);

	}
	if (tm == 1)
	{
		// 拘束
		for (i = 0; i < Total_Constraint; i++)
		{
			fscanf(fp, "%d %d %lf",
				   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0],
				   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 1],
				   &info->Value_of_Constraint[i + info->Total_Constraint_to_mesh[tm]]);
			info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] = info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];
		}
		if (Total_Constraint > 0)
		fgets(s, 256, fp);

		// 荷重
		for (i = 0; i < Total_Load; i++)
		{
			fscanf(fp, "%d %d %lf",
				   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
				   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
				   &info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
			info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] = info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];
			
			printf("Load_Node_Dir[%d]= %d Load_Node_Dir[%d]= %d Value_of_Load[%d] = %le\n",
				   (i + info->Total_Load_to_mesh[tm]) * 2 + 0, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
				   (i + info->Total_Load_to_mesh[tm]) * 2 + 1, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
				   i + info->Total_Load_to_mesh[tm], info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
		}
		if (Total_Load > 0)
		fgets(s, 256, fp);
	
		if (info->DIMENSION == 2)
		{
			int type_load, iPatch, iCoord;
			double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];
		
			for (i = 0; i < Total_DistributeForce; i++)
			{
				fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
				printf("Distibuted load number: %d\n", i);
				printf("type_load: %d iPatch: %d iCoord: %d val_Coord: %le Range_Coord: %le %le\nCoef_Dist_Load: %le %le %le\n",
						type_load, iPatch, iCoord, val_Coord, Range_Coord[0], Range_Coord[1], Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2]);
				
				// for S-IGA
				info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
				info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
				info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
				info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0] = Range_Coord[0];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1] = Range_Coord[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0] = Coeff_Dist_Load[0];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1] = Coeff_Dist_Load[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2] = Coeff_Dist_Load[2];
			}
		}
		else if (info->DIMENSION == 3)
		{
			int type_load, iPatch, iCoord, jCoord;
			double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];
		
			for (i = 0; i < Total_DistributeForce; i++)
			{
				fscanf(fp, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &jCoord, &val_Coord, &iRange_Coord[0], &iRange_Coord[1], &jRange_Coord[0], &jRange_Coord[1], &iCoeff_Dist_Load[0], &iCoeff_Dist_Load[1], &iCoeff_Dist_Load[2], &jCoeff_Dist_Load[0], &jCoeff_Dist_Load[1], &jCoeff_Dist_Load[2]);
				printf("Distibuted load number: %d\n", i);
				printf("type_load: %d iPatch: %d iCoord: %d jCoord: %d val_Coord: %le\nRange_Coord: %le %le %le %le\nCoef_Dist_Load: %le %le %le %le %le %le\n",
						type_load, iPatch, iCoord, jCoord, val_Coord,
						iRange_Coord[0], iRange_Coord[1], jRange_Coord[0], jRange_Coord[1],
						iCoeff_Dist_Load[0], iCoeff_Dist_Load[1], iCoeff_Dist_Load[2], jCoeff_Dist_Load[0], jCoeff_Dist_Load[1], jCoeff_Dist_Load[2]);
				
				// for S-IGA
				info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
				info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
				info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
				info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = jCoord;
				info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0] = iRange_Coord[0];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1] = iRange_Coord[1];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2] = jRange_Coord[0];
				info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3] = jRange_Coord[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0] = iCoeff_Dist_Load[0];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1] = iCoeff_Dist_Load[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2] = iCoeff_Dist_Load[2];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3] = jCoeff_Dist_Load[0];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4] = jCoeff_Dist_Load[1];
				info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5] = jCoeff_Dist_Load[2];
			}
		}
		if (Total_DistributeForce > 0)
		fgets(s, 256, fp);

		// ローカルパッチが属するグローバルパッチ番号(スキップ)
		fscanf(fp, "%d", &temp_i);
		fgets(s, 256, fp);

		// パッチ数(スキップ)
		fscanf(fp, "%d", &temp_i);
		fgets(s, 256, fp);
		int No_Patch = temp_i;

		// コントロールポイント数(スキップ)
		fscanf(fp, "%d", &temp_i);
		fgets(s, 256, fp);
		int Total_Control_Point = temp_i;

		// 各方向の次数
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				fscanf(fp, "%d", &temp_i);
				info->Geo_Order[i * info->DIMENSION + j] = temp_i;
				if (MAX_ORDER < temp_i)
				{
					MAX_ORDER = temp_i;
				}
				printf("Geo_Order[%d] = %d\n", i * info->DIMENSION + j, info->Geo_Order[i * info->DIMENSION + j]);
			}
		}
		fgets(s, 256, fp);

		// ノット数
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				fscanf(fp, "%d", &temp_i);
				info->Geo_No_knot[i * info->DIMENSION + j] = temp_i;
				if (MAX_KNOT < temp_i)
				{
					MAX_KNOT = temp_i;
				}
				info->Geo_Total_Knot_to_patch_dim[i * info->DIMENSION + j + 1] = info->Geo_Total_Knot_to_patch_dim[i * info->DIMENSION + j] + temp_i;
				printf("Geo_No_knot[%d] = %d\n", i * info->DIMENSION + j, info->Geo_No_knot[i * info->DIMENSION + j]);
			}
		}
		fgets(s, 256, fp);

		// 各パッチ各方向のコントロールポイント数
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				fscanf(fp, "%d", &temp_i);
				info->Geo_No_Control_point[i * info->DIMENSION + j] = temp_i;
				printf("Geo_No_Control_point[%d] = %d\n", i * info->DIMENSION + j, info->Geo_No_Control_point[i * info->DIMENSION + j]);
			}
		}
		fgets(s, 256, fp);

		for (i = 0; i < No_Patch; i++)
		{
			info->Geo_No_Control_point_in_patch[i] = 1;
		}

		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->Geo_No_Control_point_in_patch[i] *= info->Geo_No_Control_point[i * info->DIMENSION + j];
			}
		}

		for (i = 0; i < No_Patch; i++)
		{
			info->Geo_Total_Control_Point_to_patch[i + 1] = info->Geo_Total_Control_Point_to_patch[i] + info->Geo_No_Control_point_in_patch[i];
		}

		// 各パッチのコントロールポイント数
		for (i = 0; i < No_Patch; i++)
		{
			printf("Geo_No_Control_point_in_patch[%d] = %d\t", i, info->Geo_No_Control_point_in_patch[i]);
		}
		printf("\n");
	
		// パッチコネクティビティ
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->Geo_No_Control_point_in_patch[i]; j++)
			{
				fscanf(fp, "%d", &temp_i);
				info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[i] + j] = temp_i;
			}
		}
		fgets(s, 256, fp);
		
		// 境界条件(スキップ)
		int Total_Constraint, Total_Load, Total_DistributeForce;
		fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
		fgets(s, 256, fp);

		// ノットベクトルの読み込み
		for (i = 0; i < No_Patch; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				for (k = 0; k < info->Geo_No_knot[i * info->DIMENSION + j]; k++)
				{
					fscanf(fp, "%lf", &temp_d);
					info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[i * info->DIMENSION + j] + k] = temp_d;
				}
			}
		}
		fgets(s, 256, fp);
		
		int Total_Element = 0;
		for (i = 0; i < No_Patch; i++)
		{
			if (info->DIMENSION == 2)
			{
				Total_Element += (info->Geo_No_Control_point[i * info->DIMENSION + 0] - info->Geo_Order[i * info->DIMENSION + 0])
								* (info->Geo_No_Control_point[i * info->DIMENSION + 1] - info->Geo_Order[i * info->DIMENSION + 1]);
				info->Geo_No_Control_point_ON_ELEMENT[i] = (info->Geo_Order[i * info->DIMENSION + 0] + 1) * (info->Geo_Order[i * info->DIMENSION + 1] + 1);
			}
			else if (info->DIMENSION == 3)
			{
				Total_Element += (info->Geo_No_Control_point[i * info->DIMENSION + 0] - info->Geo_Order[i * info->DIMENSION + 0])
							   * (info->Geo_No_Control_point[i * info->DIMENSION + 1] - info->Geo_Order[i * info->DIMENSION + 1])
							   * (info->Geo_No_Control_point[i * info->DIMENSION + 2] - info->Geo_Order[i * info->DIMENSION + 2]);
				info->Geo_No_Control_point_ON_ELEMENT[i] = (info->Geo_Order[i * info->DIMENSION + 0] + 1) * (info->Geo_Order[i * info->DIMENSION + 1] + 1) * (info->Geo_Order[i * info->DIMENSION + 2] + 1);
			}
		}
		printf("Geo_Total_Element = %d\n", Total_Element);
		info->Geo_Total_Element_on_mesh = Total_Element;
		printf("Geo_Total_Element_on_mesh = %d\n", info->Geo_Total_Element_on_mesh);

		for (i = 0; i < No_Patch; i++)
		{
			printf("Geo_No_Control_point_ON_ELEMENT[%d] = %d\n", i, info->Geo_No_Control_point_ON_ELEMENT[i]);
		}
		
		
		// 節点座標
		for (i = 0; i < Total_Control_Point; i++)
		{
		    fscanf(fp, "%d", &temp_i);
		    for (j = 0; j < info->DIMENSION + 1; j++)
		    {
				fscanf(fp, "%lf", &temp_d);
		        double coordinate_value = temp_d;
				if (j < info->DIMENSION)
		            coordinate_value = temp_d;
	
				if (j == info->DIMENSION) 
		            coordinate_value = temp_d;
		        info->Geo_Node_Coordinate[temp_i * (info->DIMENSION + 1) + j] = coordinate_value;
		    }
		}

		for (i = 0; i < No_Patch; i++)
		{
			int cp_to_patch = info->Geo_Total_Control_Point_to_patch[i];
			for (j = 0; j < info->Geo_No_Control_point_in_patch[i]; j++)
			{
				int connectivity = info->Geo_Patch_Control_point[cp_to_patch + j] * (info->DIMENSION + 1);
				if (info->DIMENSION == 2)
				{
					info->Control_Coord_x_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity];
					info->Control_Coord_y_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity + 1];
					info->Control_Weight_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity + 2];
				}
				else if (info->DIMENSION == 3)
				{
					info->Control_Coord_x_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity];
					info->Control_Coord_y_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity + 1];
					info->Control_Coord_z_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity + 2];
					info->Control_Weight_geo[cp_to_patch + j] = info->Geo_Node_Coordinate[connectivity + 3];
				}
			}
		}
		fgets(s, 256, fp);
	}
	fclose(fp);
}


void Get_Input_3(const char *filename, information *info)
{
	// 入力ファイルから情報を読み取る関数
	char s[256];
	int temp_i, temp_i2;
	double temp_d, temp_d2, temp_d3;
	
	int i,j;

	// ローカルIGA解析モデルの要素数の取得
	int local_ele_num = info->Geo_Total_Element_on_mesh;
	
	// 要素頂点の数
	int num_vertex;
	if (info->DIMENSION == 2)
	{
		num_vertex = 9; // 2Dの場合、要素は9つの頂点を持つ
	}
	else if (info->DIMENSION == 3)
	{
		num_vertex = 27; // 3Dの場合、要素は27つの頂点を持つ
	}

	FILE *fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("warning: Could not open file %s\n", filename);
		return;
	}

	// 一行目スキップ
	fgets(s, 256, fp);

	// 要素番号, 頂点番号, ξ, η, ζの順に入力ファイルに書いてあるINPUT
	for (i = 0; i < local_ele_num; i++)
		for (j = 0; j < num_vertex; j++)
		{
			fscanf(fp, "%d %d %lf %lf %lf", &temp_i, &temp_i2, &temp_d, &temp_d2, &temp_d3);
			info->Target_Physical_Coord[(i * num_vertex + j) * info->DIMENSION + 0] = temp_d;
			info->Target_Physical_Coord[(i * num_vertex + j) * info->DIMENSION + 1] = temp_d2;
			info->Target_Physical_Coord[(i * num_vertex + j) * info->DIMENSION + 2] = temp_d3;
		}
	fclose(fp);

	// debug
	// for (i = 0; i < local_ele_num; i++)
	// {
	// 	printf("Element %d:\n", i);
	// 	for (j = 0; j < num_vertex; j++)
	// 	{
	// 		printf("  Vertex %d: (%le, %le, %le)\n", j,
	// 			   info->Target_Physical_Coord[(i * num_vertex + j) * info->DIMENSION + 0],
	// 			   info->Target_Physical_Coord[(i * num_vertex + j) * info->DIMENSION + 1],
	// 			   info->Target_Physical_Coord[(i * num_vertex + j) * info->DIMENSION + 2]);
	// 	}
	// }
}
	

// INC 等の作成
void Make_INC(information *info)
{
	int tm;

	// info->INC の計算(節点番号をξ, ηの番号で表す為の配列)
	for (tm = 0; tm < Total_mesh; tm++)
	{
		int b, B, e, h, i, j, k, l, n, o, p, q, x, y, z, ii, jj, kk, kkk, iiloc, jjloc, kkloc, r = 0;
		int type_load, iPatch, iCoord, jCoord;
		double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];
		int No_Patch = info->Total_Patch_on_mesh[tm];
		int Total_Patch_to_Now = info->Total_Patch_to_mesh[tm];
		int Total_Element = info->Total_Element_on_mesh[tm];
		int Total_Element_to_Now = info->Total_Element_to_mesh[tm];
		int Total_DistributeForce = info->Total_DistributeForce_to_mesh[tm + 1] - info->Total_DistributeForce_to_mesh[tm];

		if (info->DIMENSION == 2)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)
				{
					for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)
					{
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;

						if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1])
						{
							for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)
								{
									B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
									b = jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
									info->Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
								}
							}
							info->Element_patch[e + Total_Element_to_Now] = patch;
							info->Element_mesh[e + Total_Element_to_Now] = tm;
							e++;
						}
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (kk = 0; kk < info->No_Control_point[patch * info->DIMENSION + 2]; kk++)
				{
					for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)
					{
						for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)
						{
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 2] = kk;

							if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1] && kk >= info->Order[patch * info->DIMENSION + 2])
							{
								for (kkloc = 0; kkloc <= info->Order[patch * info->DIMENSION + 2]; kkloc++)
								{
									for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)
									{
										for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)
										{
											B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - kkloc * info->No_Control_point[patch * info->DIMENSION + 0] * info->No_Control_point[patch * info->DIMENSION + 1] - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
											b = kkloc * (info->Order[patch * info->DIMENSION + 0] + 1) * (info->Order[patch * info->DIMENSION + 1] + 1) + jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
											info->Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
										}
									}
								}
								info->Element_patch[e + Total_Element_to_Now] = patch;
								info->Element_mesh[e + Total_Element_to_Now] = tm;
								e++;
							}
							i++;
						}
					}
				}
			}
		}

		// for SS-IGA line_No_real_elementの初期化
		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->line_No_real_element[patch * info->DIMENSION + j] = 0;
			}
		}

		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->line_No_Total_element[patch * info->DIMENSION + j] = info->No_knot[patch * info->DIMENSION + j] - 2 * info->Order[patch * info->DIMENSION + j] - 1;

				for (kkk = info->Order[patch * info->DIMENSION + j]; kkk < info->No_knot[patch * info->DIMENSION + j] - info->Order[patch * info->DIMENSION + j] - 1; kkk++)
				{
					info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Order[patch * info->DIMENSION + j]]
						= info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk];
					if (info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Order[patch * info->DIMENSION + j]] != 0)
					{
						info->line_No_real_element[patch * info->DIMENSION + j]++;
					}
				}
			}
		}

		// 要素に行番号, 列番号をつける
		for (h = 0; h < Total_Element; h++)
		{
			info->Total_element_all_ID[h] = 0;
		}

		if (info->DIMENSION == 2)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (y = 0; y < info->line_No_Total_element[patch * info->DIMENSION + 1]; y++)
				{
					for (x = 0; x < info->line_No_Total_element[patch * info->DIMENSION + 0]; x++)
					{
						info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 0] = x;
						info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 1] = y;
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (z = 0; z < info->line_No_Total_element[patch * info->DIMENSION + 2]; z++)
				{
					for (y = 0; y < info->line_No_Total_element[patch * info->DIMENSION + 1]; y++)
					{
						for (x = 0; x < info->line_No_Total_element[patch * info->DIMENSION + 0]; x++)
						{
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 0] = x;
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 1] = y;
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 2] = z;
							i++;
						}
					}
				}
			}
		}

		// 必要な要素の行と列の番号を求める
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (l = 0; l < No_Patch; l++)
			{
				e = 0;
				int patch = l + Total_Patch_to_Now;
				for (k = 0; k < info->line_No_Total_element[patch * info->DIMENSION + j]; k++)
				{
					if (info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + k] != 0)
					{
						info->real_element_line[patch * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + e * info->DIMENSION + j] = k;
						e++;
					}
				}
			}
		}

		// 必要な要素列上の要素のIDを1にする
		if (info->DIMENSION == 2)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (p = 0; p < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; p++)
				{
					if (info->ENC[ele * info->DIMENSION + 0] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + p * info->DIMENSION + 0])
					{
						for (q = 0; q < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; q++)
						{
							if (info->ENC[ele * info->DIMENSION + 1] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + q * info->DIMENSION + 1])
							{
								info->Total_element_all_ID[n]++;
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Total_element_all_ID[n] == 1)
				{
					info->real_element[r + info->real_Total_Element_to_mesh[tm]] = ele;
					info->real_El_No_on_mesh[tm * info->Total_Element_to_mesh[Total_mesh] + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->line_No_real_element[patch * info->DIMENSION + 0] * info->line_No_real_element[patch * info->DIMENSION + 1];
			}
			info->real_Total_Element_on_mesh[tm] = real_Total_Element;
			info->real_Total_Element_to_mesh[tm + 1] = info->real_Total_Element_to_mesh[tm] + real_Total_Element;
		}
		else if (info->DIMENSION == 3)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (o = 0; o < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; o++)
				{
					if (info->ENC[ele * info->DIMENSION + 0] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + o * info->DIMENSION + 0])
					{
						for (p = 0; p < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; p++)
						{
							if (info->ENC[ele * info->DIMENSION + 1] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + p * info->DIMENSION + 1])
							{
								for (q = 0; q < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 2]; q++)
								{
									if (info->ENC[ele * info->DIMENSION + 2] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + q * info->DIMENSION + 2])
									{
										info->Total_element_all_ID[n]++;
									}
								}
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Total_element_all_ID[n] == 1)
				{
					info->real_element[r + info->real_Total_Element_to_mesh[tm]] = ele;
					info->real_El_No_on_mesh[tm * info->Total_Element_to_mesh[Total_mesh] + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->line_No_real_element[patch * info->DIMENSION + 0] * info->line_No_real_element[patch * info->DIMENSION + 1] * info->line_No_real_element[patch * info->DIMENSION + 2];
			}
			info->real_Total_Element_on_mesh[tm] = real_Total_Element;
			info->real_Total_Element_to_mesh[tm + 1] = info->real_Total_Element_to_mesh[tm] + real_Total_Element;
		}

		// 形状出力のみ
		if (info->c.GEOMETRY_ONLY_OUTPUT == 1)
			continue;

		// 分布荷重
		for (i = 0; i < Total_DistributeForce; i++)
		{
			if (info->DIMENSION == 2)
			{
				type_load = info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iPatch = info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iCoord = info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				val_Coord = info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0];
				iRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1];
				iCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0];
				iCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1];
				iCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2];

				setDistLoad(tm, iPatch, iCoord, val_Coord, iRange_Coord, type_load, iCoeff_Dist_Load, info);
				// setDistLoad_infinite_plate_with_hole(tm, iPatch, iCoord, val_Coord, iRange_Coord, type_load, iCoeff_Dist_Load, info);
			}
			else if (info->DIMENSION == 3)
			{
				type_load = info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iPatch = info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iCoord = info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				jCoord = info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				val_Coord = info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0];
				iRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1];
				jRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2];
				jRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3];
				iCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0];
				iCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1];
				iCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2];
				jCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3];
				jCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4];
				jCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5];

				setDistLoad(tm, iPatch, iCoord, jCoord, val_Coord, iRange_Coord, jRange_Coord, type_load, iCoeff_Dist_Load, jCoeff_Dist_Load, info);
			}
		}
	}
}


// SSIGAローカル形状表現のためのINC作成
void Make_INC_for_local_geometric(information *info)
{
	int tm;
	int Total_mesh_geo = 1;		// ジオメトリ表現のメッシュは一個だけ

	for (tm = 0; tm < Total_mesh_geo; tm++)
	{
		int b, B, h, n, o, p, q, x, y, z;
		int e, i, j, k, l, ii, jj, kk, kkk, iiloc, jjloc, kkloc, r = 0;
		int No_Patch = info->Geo_Total_patch_on_mesh;
		int Total_Patch_to_Now = info->Geo_Total_patch_to_mesh;
		int Total_Element = info->Geo_Total_Element_on_mesh;
		int Total_Element_to_Now = info->Geo_Total_Element_to_mesh;


		if (info->DIMENSION == 2)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (jj = 0; jj < info->Geo_No_Control_point[patch * info->DIMENSION + 1]; jj++)
				{
					for (ii = 0; ii < info->Geo_No_Control_point[patch * info->DIMENSION + 0]; ii++)
					{
						info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
						info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;

						if (ii >= info->Geo_Order[patch * info->DIMENSION + 0] && jj >= info->Geo_Order[patch * info->DIMENSION + 1])
						{
							for (jjloc = 0; jjloc <= info->Geo_Order[patch * info->DIMENSION + 1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= info->Geo_Order[patch * info->DIMENSION + 0]; iiloc++)
								{
									B = info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i - jjloc * info->Geo_No_Control_point[patch * info->DIMENSION + 0] - iiloc];
									b = jjloc * (info->Geo_Order[patch * info->DIMENSION + 0] + 1) + iiloc;
									info->Geo_Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
								}
							}
							info->Geo_Element_patch[e + Total_Element_to_Now] = patch;
							info->Geo_Element_mesh[e + Total_Element_to_Now] = tm;
							e++;
						}
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (kk = 0; kk < info->Geo_No_Control_point[patch * info->DIMENSION + 2]; kk++)
				{
					for (jj = 0; jj < info->Geo_No_Control_point[patch * info->DIMENSION + 1]; jj++)
					{
						for (ii = 0; ii < info->Geo_No_Control_point[patch * info->DIMENSION + 0]; ii++)
						{
							info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
							info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;
							info->Geo_INC[patch * info->Geo_Total_Control_Point_on_mesh * info->DIMENSION + info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 2] = kk;

							if (ii >= info->Geo_Order[patch * info->DIMENSION + 0] && jj >= info->Geo_Order[patch * info->DIMENSION + 1] && kk >= info->Geo_Order[patch * info->DIMENSION + 2])
							{
								for (kkloc = 0; kkloc <= info->Geo_Order[patch * info->DIMENSION + 2]; kkloc++)
								{
									for (jjloc = 0; jjloc <= info->Geo_Order[patch * info->DIMENSION + 1]; jjloc++)
									{
										for (iiloc = 0; iiloc <= info->Geo_Order[patch * info->DIMENSION + 0]; iiloc++)
										{
											B = info->Geo_Patch_Control_point[info->Geo_Total_Control_Point_to_patch[patch] + i - kkloc * info->Geo_No_Control_point[patch * info->DIMENSION + 0] * info->Geo_No_Control_point[patch * info->DIMENSION + 1] - jjloc * info->Geo_No_Control_point[patch * info->DIMENSION + 0] - iiloc];
											b = kkloc * (info->Geo_Order[patch * info->DIMENSION + 0] + 1) * (info->Geo_Order[patch * info->DIMENSION + 1] + 1) + jjloc * (info->Geo_Order[patch * info->DIMENSION + 0] + 1) + iiloc;
											info->Geo_Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
										}
									}
								}
								info->Geo_Element_patch[e + Total_Element_to_Now] = patch;
								info->Geo_Element_mesh[e + Total_Element_to_Now] = tm;
								e++;
							}
							i++;
						}
					}
				}
			}
		}

		// for SS-IGA line_No_real_elementの初期化
		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->Geo_line_No_real_element[patch * info->DIMENSION + j] = 0;
			}
		}

		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->Geo_line_No_Total_element[patch * info->DIMENSION + j] = info->Geo_No_knot[patch * info->DIMENSION + j] - 2 * info->Geo_Order[patch * info->DIMENSION + j] - 1;

				for (kkk = info->Geo_Order[patch * info->DIMENSION + j]; kkk < info->Geo_No_knot[patch * info->DIMENSION + j] - info->Geo_Order[patch * info->DIMENSION + j] - 1; kkk++)
				{
					info->Geo_difference[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Geo_Order[patch * info->DIMENSION + j]]
						= info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk + 1] - info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk];
					if (info->Geo_difference[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Geo_Order[patch * info->DIMENSION + j]] != 0)
					{
						info->Geo_line_No_real_element[patch * info->DIMENSION + j]++;
					}
				}
			}
		}

		// 要素に行番号, 列番号をつける
		for (h = 0; h < Total_Element; h++)
		{
			info->Geo_Total_element_all_ID[h] = 0;
		}

		if (info->DIMENSION == 2)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (y = 0; y < info->Geo_line_No_Total_element[patch * info->DIMENSION + 1]; y++)
				{
					for (x = 0; x < info->Geo_line_No_Total_element[patch * info->DIMENSION + 0]; x++)
					{
						info->Geo_ENC[(i + info->Geo_Total_Element_to_mesh) * info->DIMENSION + 0] = x;
						info->Geo_ENC[(i + info->Geo_Total_Element_to_mesh) * info->DIMENSION + 1] = y;
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (z = 0; z < info->Geo_line_No_Total_element[patch * info->DIMENSION + 2]; z++)
				{
					for (y = 0; y < info->Geo_line_No_Total_element[patch * info->DIMENSION + 1]; y++)
					{
						for (x = 0; x < info->Geo_line_No_Total_element[patch * info->DIMENSION + 0]; x++)
						{
							info->Geo_ENC[(i + info->Geo_Total_Element_to_mesh) * info->DIMENSION + 0] = x;
							info->Geo_ENC[(i + info->Geo_Total_Element_to_mesh) * info->DIMENSION + 1] = y;
							info->Geo_ENC[(i + info->Geo_Total_Element_to_mesh) * info->DIMENSION + 2] = z;
							i++;
						}
					}
				}
			}
		}

		// 必要な要素の行と列の番号を求める
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (l = 0; l < No_Patch; l++)
			{
				e = 0;
				int patch = l + Total_Patch_to_Now;
				for (k = 0; k < info->Geo_line_No_Total_element[patch * info->DIMENSION + j]; k++)
				{
					if (info->Geo_difference[info->Geo_Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + k] != 0)
					{
						info->Geo_real_element_line[patch * (info->Geo_Total_Element_to_mesh * info->DIMENSION) + e * info->DIMENSION + j] = k;
						e++;
					}
				}
			}
		}

		// 必要な要素列上の要素のIDを1にする
		if (info->DIMENSION == 2)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (p = 0; p < info->Geo_line_No_real_element[info->Geo_Element_patch[ele] * info->DIMENSION + 0]; p++)
				{
					if (info->Geo_ENC[ele * info->DIMENSION + 0] == info->Geo_real_element_line[info->Geo_Element_patch[ele] * (info->Geo_Total_Element_to_mesh * info->DIMENSION) + p * info->DIMENSION + 0])
					{
						for (q = 0; q < info->Geo_line_No_real_element[info->Geo_Element_patch[ele] * info->DIMENSION + 1]; q++)
						{
							if (info->Geo_ENC[ele * info->DIMENSION + 1] == info->Geo_real_element_line[info->Geo_Element_patch[ele] * (info->Geo_Total_Element_to_mesh * info->DIMENSION) + q * info->DIMENSION + 1])
							{
								info->Geo_Total_element_all_ID[n]++;
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Geo_Total_element_all_ID[n] == 1)
				{
					info->Geo_real_element[r + info->Geo_real_Total_Element_to_mesh] = ele;
					info->Geo_real_El_No_on_mesh[tm * info->Geo_Total_Element_to_mesh + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->Geo_line_No_real_element[patch * info->DIMENSION + 0] * info->Geo_line_No_real_element[patch * info->DIMENSION + 1];
			}
			info->Geo_real_Total_Element_on_mesh = real_Total_Element;
			info->Geo_real_Total_Element_to_mesh = info->Geo_real_Total_Element_to_mesh + real_Total_Element;
		}
		else if (info->DIMENSION == 3)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (o = 0; o < info->Geo_line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; o++)
				{
					if (info->Geo_ENC[ele * info->DIMENSION + 0] == info->Geo_real_element_line[info->Geo_Element_patch[ele] * (info->Geo_Total_Element_to_mesh * info->DIMENSION) + o * info->DIMENSION + 0])
					{
						for (p = 0; p < info->Geo_line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; p++)
						{
							if (info->Geo_ENC[ele * info->DIMENSION + 1] == info->Geo_real_element_line[info->Geo_Element_patch[ele] * (info->Geo_Total_Element_to_mesh * info->DIMENSION) + p * info->DIMENSION + 1])
							{
								for (q = 0; q < info->Geo_line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 2]; q++)
								{
									if (info->Geo_ENC[ele * info->DIMENSION + 2] == info->Geo_real_element_line[info->Geo_Element_patch[ele] * (info->Geo_Total_Element_to_mesh * info->DIMENSION) + q * info->DIMENSION + 2])
									{
										info->Geo_Total_element_all_ID[n]++;
									}
								}
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Geo_Total_element_all_ID[n] == 1)
				{
					info->Geo_real_element[r + info->Geo_real_Total_Element_to_mesh] = ele;
					info->Geo_real_El_No_on_mesh[tm * info->Geo_Total_Element_to_mesh + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->Geo_line_No_real_element[patch * info->DIMENSION + 0] * info->Geo_line_No_real_element[patch * info->DIMENSION + 1] * info->Geo_line_No_real_element[patch * info->DIMENSION + 2];
			}
			info->Geo_real_Total_Element_on_mesh = real_Total_Element;
			info->Geo_real_Total_Element_to_mesh = info->Geo_real_Total_Element_to_mesh + real_Total_Element;
		}

	}
}



void setDistLoad(int current_mesh, int patch, int coord, double target_knot, double *range, int type_load, double *dist_load_coeff, information *info)
{
	Make_gauss_array(info);

	// check coordinates and set load surface
	int target_coord = -1;
	if (coord == 0)
		target_coord = 1;
	else if (coord == 1)
		target_coord = 0;
	else
	{
		printf("Error, incorrect input data at distributed load.\n");
		exit(1);
	}

	// set target ENC
	int current_patch = patch + info->Total_Patch_to_mesh[current_mesh];
	int id = current_patch * info->DIMENSION + target_coord;
	int start_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->Order[id]];
	int end_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->No_knot[id] - info->Order[id] - 1];

	bool isStart = false;
	int target_ENC = -1;
	if (fabs(target_knot - start_knot) < MERGE_ERROR)
	{
		target_ENC = 0;
		isStart = true;
	}
	else if (fabs(target_knot - end_knot) < MERGE_ERROR)
	{
		target_ENC = info->No_Control_point[id] - info->Order[id] - 1;
		isStart = false;
	}
	else
	{
		printf("Error, incorrect input data at distributed load.\ntarget_knot is invalid.\n");
		exit(1);
	}
	
	// set coord ENC
	vector<int> ENC(2);
	bool errorFlag[2] = {false, false};
	int id_coord = current_patch * info->DIMENSION + coord;
	int knot_id_coord = info->Total_Knot_to_patch_dim[id_coord] + info->Order[id_coord];
	for (int i = 0; i < info->No_Control_point[id_coord] - info->Order[id_coord]; i++)
	{
		if (fabs(range[0] - info->Position_Knots[knot_id_coord + i]) < MERGE_ERROR)
		{
			ENC[0] = i;
			errorFlag[0] |= true;
		}
		if (fabs(range[1] - info->Position_Knots[knot_id_coord + i + 1]) < MERGE_ERROR)
		{
			ENC[1] = i;
			errorFlag[1] |= true;
		}
	}
	if (!errorFlag[0] && !errorFlag[1])
	{
		printf("Error, incorrect input data at distributed load.\nChange the range of the distributed load.\n");
		exit(1);
	}

	// serach element
	vector<int> ele_list;
	for (int i = 0; i < info->Total_Element_on_mesh[current_mesh]; i++)
	{
		int e = i + info->Total_Element_to_mesh[current_mesh];
		if (info->Element_patch[e] == current_patch)
		{
			bool isTargetEle = true;
			for (int j = 0; j < info->DIMENSION; j++)
			{
				int current_ENC = info->ENC[e * info->DIMENSION + j];

				// target coord
				if (j == target_coord)
				{
					if (current_ENC == target_ENC)
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}

				// check range
				else if (j == coord)
				{
					if (ENC[0] <= current_ENC && current_ENC <= ENC[1])
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}
			}

			if (isTargetEle)
				ele_list.emplace_back(e);
		}
	}

	// parametric coordinates
	int gp_1d = info->c.NUM_GAUSS_POINTS;
	vector<double> gp_para(gp_1d * info->DIMENSION);
	vector<double> gp_w(gp_1d);
	for (int i = 0; i < gp_1d; i++)
	{
		gp_para[i * info->DIMENSION + coord] = info->gauss_point_1D[i];
		gp_para[i * info->DIMENSION + target_coord] = isStart ? -1.0 : 1.0;
		gp_w[i] = info->gauss_w_1D[i];
	}

	// integration for distributed load
	for (size_t i = 0; i < ele_list.size(); i++)
	{
		int e = ele_list[i];
		for (int j = 0; j < gp_1d; j++)
		{
			double *para = gp_para.data() + j * info->DIMENSION;
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
			shape_and_dshape(R.data(), dR.data(), para, e, true, info);

			// jacobian matrix
			vector<double> jac(info->DIMENSION, 0.0);
			for (int k = 0; k < info->DIMENSION; k++)
				for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; m++)
					for (int l = 0; l < info->DIMENSION; l++)
						if (l == coord)
							jac[k] += dR[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];

			// jacobian ||normal_J||
			vector<double> normal_J(info->DIMENSION);
			normal_J[0] = jac[1];
			normal_J[1] = -jac[0];
			double J = 0.0;
			for (int k = 0; k < info->DIMENSION; k++)
			{
				J += normal_J[k] * normal_J[k];
			}
			J = sqrt(J);

			// integrand
			double *c = dist_load_coeff;
			vector<double> current_knot(info->DIMENSION);
			ShapeFunc_from_paren(current_knot.data(), para, coord, e, info);
			double x = current_knot[coord];

			double val = (c[0] + c[1] * x + c[2] * x * x);

			// type 0, 1
			if (type_load < 3)
			{
				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + type_load;
					info->Equivalent_Nodal_Force[index] += coeff * R[k];
				}
			}

			// type 3 (normal direction)
			else if (type_load == 2)
			{
				vector<double> normal(info->DIMENSION);
				normal[0] = normal_J[0] / J;
				normal[1] = normal_J[1] / J;

				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION;
					info->Equivalent_Nodal_Force[index + 0] += coeff * normal[0] * R[k];
					info->Equivalent_Nodal_Force[index + 1] += coeff * normal[1] * R[k];
				}
			}

			else
			{
				printf("Error, incorrect input data at distributed load.\nChange the type of the distributed load.\n");
				exit(1);
			}
		}
	}
}


void setDistLoad(int current_mesh, int patch, int coord_i, int coord_j, double target_knot, double *range_i, double *range_j, int type_load, double *dist_load_coeff_i, double *dist_load_coeff_j, information *info)
{
	Make_gauss_array(info);

	// check coordinates and set load surface
	int target_coord = -1;
	if (coord_i == 0 && coord_j == 1)
		target_coord = 2;
	else if (coord_i == 1 && coord_j == 2)
		target_coord = 0;
	else if (coord_i == 2 && coord_j == 0)
		target_coord = 1;
	else
	{
		printf("Error, incorrect input data at distributed load.\n");
		exit(1);
	}

	// set target ENC
	int current_patch = patch + info->Total_Patch_to_mesh[current_mesh];
	int id = current_patch * info->DIMENSION + target_coord;
	int start_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->Order[id]];
	int end_knot = info->Position_Knots[info->Total_Knot_to_patch_dim[id] + info->No_knot[id] - info->Order[id] - 1];

	// check local or global flag
	bool isLocal = false;
	if (current_patch >= info->Total_Patch_to_mesh[1])
		isLocal = true;

	bool isStart = false;
	int target_ENC = -1;
	if (fabs(target_knot - start_knot) < MERGE_ERROR)
	{
		target_ENC = 0;
		isStart = true;
	}
	else if (fabs(target_knot - end_knot) < MERGE_ERROR)
	{
		target_ENC = info->No_Control_point[id] - info->Order[id] - 1;
		isStart = false;
	}
	else
	{
		printf("Error, incorrect input data at distributed load.\ntarget_knot is invalid.\n");
		exit(1);
	}

	// set coord_i ENC
	vector<int> i_ENC(2);
	bool errorFlag_i[2] = {false, false};
	int id_i = current_patch * info->DIMENSION + coord_i;
	int knot_id_i = info->Total_Knot_to_patch_dim[id_i] + info->Order[id_i];
	for (int i = 0; i < info->No_Control_point[id_i] - info->Order[id_i]; i++)
	{
		if (fabs(range_i[0] - info->Position_Knots[knot_id_i + i]) < MERGE_ERROR)
		{
			i_ENC[0] = i;
			errorFlag_i[0] |= true;
		}
		if (fabs(range_i[1] - info->Position_Knots[knot_id_i + i + 1]) < MERGE_ERROR)
		{
			i_ENC[1] = i;
			errorFlag_i[1] |= true;
		}
	}
	if (!errorFlag_i[0] && !errorFlag_i[1])
	{
		printf("Error, incorrect input data at distributed load.\nChange the range of the distributed load.\n");
		exit(1);
	}

	// set coord_j ENC
	vector<int> j_ENC(2);
	bool errorFlag_j[2] = {false, false};
	int id_j = current_patch * info->DIMENSION + coord_j;
	int knot_id_j = info->Total_Knot_to_patch_dim[id_j] + info->Order[id_j];
	for (int i = 0; i < info->No_Control_point[id_j] - info->Order[id_j]; i++)
	{
		if (fabs(range_j[0] - info->Position_Knots[knot_id_j + i]) < MERGE_ERROR)
		{
			j_ENC[0] = i;
			errorFlag_j[0] |= true;
		}
		if (fabs(range_j[1] - info->Position_Knots[knot_id_j + i + 1]) < MERGE_ERROR)
		{
			j_ENC[1] = i;
			errorFlag_j[1] |= true;
		}
	}
	if (!errorFlag_j[0] && !errorFlag_j[1])
	{
		printf("Error, incorrect input data at distributed load.\nChange the range of the distributed load.\n");
		exit(1);
	}

	// serach element
	vector<int> ele_list;
	for (int i = 0; i < info->Total_Element_on_mesh[current_mesh]; i++)
	{
		int e = i + info->Total_Element_to_mesh[current_mesh];
		if (info->Element_patch[e] == current_patch)
		{
			bool isTargetEle = true;
			for (int j = 0; j < info->DIMENSION; j++)
			{
				int current_ENC = info->ENC[e * info->DIMENSION + j];

				// target coord
				if (j == target_coord)
				{
					if (current_ENC == target_ENC)
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}

				// check range i
				else if (j == coord_i)
				{
					if (i_ENC[0] <= current_ENC && current_ENC <= i_ENC[1])
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}

				// check range j
				else if (j == coord_j)
				{
					if (j_ENC[0] <= current_ENC && current_ENC <= j_ENC[1])
						isTargetEle &= true;
					else
					{
						isTargetEle &= false;
						break;
					}
				}
			}

			if (isTargetEle)
				ele_list.emplace_back(e);
		}
	}

	// parametric coordinates
	int gp_1d = info->c.NUM_GAUSS_POINTS;
	int gp_2d = info->c.NUM_GAUSS_POINTS * info->c.NUM_GAUSS_POINTS;
	vector<double> gp_para(gp_2d * info->DIMENSION);
	vector<double> gp_w(gp_2d);
	for (int i = 0; i < gp_1d; i++)
		for (int j = 0; j < gp_1d; j++)
		{
			gp_para[(i * gp_1d + j) * info->DIMENSION + coord_i] = info->gauss_point_1D[i];
			gp_para[(i * gp_1d + j) * info->DIMENSION + coord_j] = info->gauss_point_1D[j];
			gp_para[(i * gp_1d + j) * info->DIMENSION + target_coord] = isStart ? -1.0 : 1.0;
			gp_w[i * gp_1d + j] = info->gauss_w_1D[i] * info->gauss_w_1D[j];
		}

	// integration for distributed load for global coordinate direction
	for (size_t i = 0; i < ele_list.size(); i++)
	{
		int e = ele_list[i];
		for (int j = 0; j < gp_2d; j++)
		{
			double *para = gp_para.data() + j * info->DIMENSION;
			vector<double> R(MAX_NO_CP_ON_ELEMENT);
			vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

			vector<double> jac_i(info->DIMENSION, 0.0);
			vector<double> jac_j(info->DIMENSION, 0.0); 
			double J = 0.0;
			vector<double> normal_J(info->DIMENSION);

			if (!isLocal)
			{
				shape_and_dshape(R.data(), dR.data(), para, e, true, info);

				// jacobian matrix
				for (int k = 0; k < info->DIMENSION; k++)
					for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; m++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							if (l == coord_i)
								jac_i[k] += dR[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];
							else if (l == coord_j)
								jac_j[k] += dR[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];
						}

				// jacobian ||cross(jac_i, jac_j)||
				for (int k = 0; k < info->DIMENSION; k++)
				{
					normal_J[k] = jac_i[(k + 1) % info->DIMENSION] * jac_j[(k + 2) % info->DIMENSION] - jac_i[(k + 2) % info->DIMENSION] * jac_j[(k + 1) % info->DIMENSION];
					J += normal_J[k] * normal_J[k];
				}
				J = sqrt(J);
			}
			else
			{
				vector<double> a(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> b(info->DIMENSION * info->DIMENSION, 0.0);
				vector<double> c(info->DIMENSION * info->DIMENSION, 0.0);

				vector<double> para_geo(info->DIMENSION); // ローカル形状要素パラメータ座標
				vector<double> para_glo(info->DIMENSION); // グローバル要素パラメータ座標
				vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
				vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
				vector<double> dR_geo(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
				vector<double> dR_glo(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

				int ele_geo = trans_local_para_to_local_geo_para(e, para, para_geo.data(), info);
				int ele_glo = trans_local_para_to_global_para(e, para, para_glo.data(), info);

				Bspline_shape_and_dshape(R.data(), dR.data(), para, e, true, info);
				geo_shape_and_dshape(R_geo.data(), dR_geo.data(), para_geo.data(), ele_geo, false, info);
				shape_and_dshape(R_glo.data(), dR_glo.data(), para_glo.data(), ele_glo, false, info);

				// 変位要素空間の微分によって生じる定数
				for (int i = 0; i < info->DIMENSION; i++)
				{
					double geo_coeff = dShapeFunc_from_paren(i, e, info);
					for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; j++)
						dR_geo[j * info->DIMENSION + i] *= geo_coeff;
				}
				
				// jacobian matrix
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; m++)
						{
							a[k * info->DIMENSION + l] += dR_geo[m * info->DIMENSION + l] * info->Geo_Node_Coordinate[info->Geo_Controlpoint_of_Element[ele_geo * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];
						}
				for (int k = 0; k < info->DIMENSION; k++)
					for (int l = 0; l < info->DIMENSION; l++)
						for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_glo]]; m++)
						{
							b[k * info->DIMENSION + l] += dR_glo[m * info->DIMENSION + l] * info->Node_Coordinate[info->Controlpoint_of_Element[ele_glo * MAX_NO_CP_ON_ELEMENT + m] * (info->DIMENSION + 1) + k];
						}

        		for (int k = 0; k < info->DIMENSION; k++)
        		{
					for (int l = 0; l < info->DIMENSION; l++)
        		    {
						for (int m = 0; m < info->DIMENSION; m++)
        		        {
        		            c[k * info->DIMENSION + l] += b[k * info->DIMENSION + m] * a[m * info->DIMENSION + l];
        		        }
        		    }
        		}

				for (int k = 0; k < info->DIMENSION; k++)
				{
					jac_i[k] = c[k * info->DIMENSION + coord_i];
					jac_j[k] = c[k * info->DIMENSION + coord_j];
				}
			
				// jacobian ||cross(jac_i, jac_j)||
				for (int k = 0; k < info->DIMENSION; k++)
				{
					normal_J[k] = jac_i[(k + 1) % info->DIMENSION] * jac_j[(k + 2) % info->DIMENSION] - jac_i[(k + 2) % info->DIMENSION] * jac_j[(k + 1) % info->DIMENSION];
					J += normal_J[k] * normal_J[k];
				}
				J = sqrt(J);
			}

			// integrand
			double *ci = dist_load_coeff_i;
			double *cj = dist_load_coeff_j;
			vector<double> current_knot(info->DIMENSION);
			ShapeFunc_from_paren(current_knot.data(), para, coord_i, e, info);
			ShapeFunc_from_paren(current_knot.data(), para, coord_j, e, info);
			double x = current_knot[coord_i];
			double y = current_knot[coord_j];

			double val = (ci[0] + ci[1] * x + ci[2] * x * x) * (cj[0] + cj[1] * y + cj[2] * y * y);

			// type 0, 1, 2
			if (type_load < 3)
			{
				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + type_load;
					info->Equivalent_Nodal_Force[index] += coeff * R[k];
				}
			}

			// type 3 (normal direction)
			else if (type_load == 3)
			{
				vector<double> normal(info->DIMENSION);
				normal[0] = normal_J[0] / J;
				normal[1] = normal_J[1] / J;
				normal[2] = normal_J[2] / J;

				double coeff = val * J * gp_w[j];
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
				{
					int index = info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION;
					info->Equivalent_Nodal_Force[index + 0] += coeff * normal[0] * R[k];
					info->Equivalent_Nodal_Force[index + 1] += coeff * normal[1] * R[k];
					info->Equivalent_Nodal_Force[index + 2] += coeff * normal[2] * R[k];
				}
			}

			else
			{
				printf("Error, incorrect input data at distributed load.\nChange the type of the distributed load.\n");
				exit(1);
			}
		}
	}

	// exit(0);
}

// 要素のパラメータ座標からパッチのパラメータ座標を計算する関数
void trans_ele_patch_coord(double* xi_patch, const double* xi_elem, int patch_num, int ele_num, information* info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int knot_start = info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i];
		double knot_front = info->Position_Knots[knot_start];
		double knot_back  = info->Position_Knots[knot_start + 1];
		xi_patch[i] = (xi_elem[i] + 1.0) / 2.0 * (knot_back - knot_front) + knot_front;
	}
}

// ローカル形状表現の要素のパラメータ座標からローカルパッチのパラメータ座標を計算する関数
void geo_trans_ele_patch_coord(double* xi_patch, const double* xi_elem, int patch_num, int ele_num, information* info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int knot_start = info->Geo_Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->Geo_ENC[ele_num * info->DIMENSION + i];
		double knot_front = info->Geo_Position_Knots[knot_start];
		double knot_back  = info->Geo_Position_Knots[knot_start + 1];
		xi_patch[i] = (xi_elem[i] + 1.0) / 2.0 * (knot_back - knot_front) + knot_front;
	}
}


// search element
int ele_check(int patch_n, double *para_coord, information *info)
{
	int line[MAX_DIMENSION] = {0};
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int line_counter = 0;
		for (int j = 0; j < info->No_Control_point[patch_n * info->DIMENSION + i] - info->Order[patch_n * info->DIMENSION + i]; j++)
		{
			double knot_front = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Order[patch_n * info->DIMENSION + i] + j];
			double knot_back = info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Order[patch_n * info->DIMENSION + i] + j + 1];

			if (fabs(knot_front - para_coord[i]) < MERGE_ERROR)
			{
				para_coord[i] = knot_front;
				line[i] = line_counter;
				break;
			}

			if (fabs(knot_back - para_coord[i]) < MERGE_ERROR)
			{
				para_coord[i] = knot_back;
				line[i] = line_counter;
				break;
			}

			if (knot_front < para_coord[i] && para_coord[i] < knot_back)
			{
				line[i] = line_counter;
				break;
			}

			line_counter++;
		}
	}

	int ele = 0;
	if (info->DIMENSION == 2)
		ele = line[0] + line[1] * info->line_No_Total_element[patch_n * info->DIMENSION + 0];
	else if (info->DIMENSION == 3)
		ele = line[0] + line[1] * info->line_No_Total_element[patch_n * info->DIMENSION + 0] + line[2] * info->line_No_Total_element[patch_n * info->DIMENSION + 0] * info->line_No_Total_element[patch_n * info->DIMENSION + 1];

	for (int i = 0; i < patch_n; i++)
	{
		int ele_in_patch = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			ele_in_patch *= (info->No_Control_point[i * info->DIMENSION + j] - info->Order[i * info->DIMENSION + j]);
		ele += ele_in_patch;
	}

	return ele;
}


// ローカルのパラメータ空間座標を与えると、ローカル形状表現の要素を得る
int geo_ele_check(int patch_n, double *para_coord, information *info)
{
	int line[MAX_DIMENSION];
	for (int i = 0; i < info->DIMENSION; i++)
	{
		int line_counter = 0;
		for (int j = 0; j < info->Geo_No_Control_point[patch_n * info->DIMENSION + i] - info->Geo_Order[patch_n * info->DIMENSION + i]; j++)
		{
			double knot_front = info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Geo_Order[patch_n * info->DIMENSION + i] + j];
			double knot_back = info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch_n * info->DIMENSION + i] + info->Geo_Order[patch_n * info->DIMENSION + i] + j + 1];

			if (j == 0 && fabs(knot_front - para_coord[i]) < MERGE_ERROR)
			{
				para_coord[i] = knot_front;
				line[i] = line_counter;
				break;
			}

			if (fabs(knot_back - para_coord[i]) < MERGE_ERROR)
			{
				para_coord[i] = knot_back;
				line[i] = line_counter;
				break;
			}

			if (knot_front < para_coord[i] && para_coord[i] < knot_back)
			{
				line[i] = line_counter;
				break;
			}

			line_counter++;
		}
	}

	int ele = 0;
	if (info->DIMENSION == 2)
		ele = line[0] + line[1] * info->Geo_line_No_Total_element[patch_n * info->DIMENSION + 0];
	else if (info->DIMENSION == 3)
		ele = line[0] + line[1] * info->Geo_line_No_Total_element[patch_n * info->DIMENSION + 0] + line[2] * info->Geo_line_No_Total_element[patch_n * info->DIMENSION + 0] * info->Geo_line_No_Total_element[patch_n * info->DIMENSION + 1];

	for (int i = 0; i < patch_n; i++)
	{
		int ele_in_patch = 1;
		for (int j = 0; j < info->DIMENSION; j++)
			ele_in_patch *= (info->Geo_No_Control_point[i * info->DIMENSION + j] - info->Geo_Order[i * info->DIMENSION + j]);
		ele += ele_in_patch;
	}

	return ele;
}


void tilde_coord(double *xi_tilde, double *xi, int patch_num, int ele_num, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		xi_tilde[i] = -1.0 + 2.0 * (xi[i] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i]]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Order[patch_num * info->DIMENSION + i] + info->ENC[ele_num * info->DIMENSION + i]]);
	}
}


// ローカル形状表現用の要素パラメータ座標変換
void geo_tilde_coord(double *xi_tilde, double *xi, int patch_num, int ele_num, information *info)
{
	for (int i = 0; i < info->DIMENSION; i++)
	{
		xi_tilde[i] = -1.0 + 2.0 * (xi[i] - info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Geo_Order[patch_num * info->DIMENSION + i] + info->Geo_ENC[ele_num * info->DIMENSION + i]]) / (info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Geo_Order[patch_num * info->DIMENSION + i] + info->Geo_ENC[ele_num * info->DIMENSION + i] + 1] - info->Geo_Position_Knots[info->Geo_Total_Knot_to_patch_dim[patch_num * info->DIMENSION + i] + info->Geo_Order[patch_num * info->DIMENSION + i] + info->Geo_ENC[ele_num * info->DIMENSION + i]]);
	}
}


// 昇順ソート
void sort(int total, int *element_n_point)
{
	int i, j;
	int temp;

	for (i = 0; i < total; i++)
	{
		for (j = i + 1; j < total; j++)
		{
			if (element_n_point[i] > element_n_point[j])
			{
				temp = element_n_point[i];
				element_n_point[i] = element_n_point[j];
				element_n_point[j] = temp;
			}
		}
	}
}


// gauss points
void Make_gauss_array(information *info)
{
	if (!info->isSetGaussArray)
	{
		Gauss_point(info, info->c.NUM_GAUSS_POINTS, info->gauss_w, info->gauss_point, info->gauss_w_1D, info->gauss_point_1D);
		if (info->c.USE_EXTENDED_QUADRATURE == 1)
			Gauss_point(info, info->c.NUM_GAUSS_POINTS_EXTENDED, info->gauss_w_ex, info->gauss_point_ex, info->gauss_w_1D_ex, info->gauss_point_1D_ex);
		info->isSetGaussArray = true;
	}
}


void Gauss_point(information *info, int temp_GP_1D, double *w_ptr, double *point_ptr, double *w_1D_ptr, double *point_1D_ptr)
{
	if (info->DIMENSION == 2)
	{
		if (temp_GP_1D == 3)
		{
			double G1 = pow((3.0 / 5.0), 0.5);
			double G_vec[3] = {-G1, 0.0, G1};
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;
			double w_vec[3] = {w2, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 4)
		{
			double A = pow((6.0 / 5.0), 0.5);
			double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
			double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
			double G_vec[4] = {-G2, -G1, G1, G2};
			double B = pow(30.0, 0.5);
			double w1 = (18.0 + B) / 36.0;
			double w2 = (18.0 - B) / 36.0;
			double w_vec[4] = {w2, w1, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 5)
		{
			double A = pow((10.0 / 7.0), 0.5);
			double G1 = pow((5.0 - 2.0 * A), 0.5) / 3.0;
			double G2 = pow((5.0 + 2.0 * A), 0.5) / 3.0;
			double G_vec[5] = {-G2, -G1, 0.0, G1, G2};
			double B = pow(70.0, 0.5);
			double w1 = 128.0 / 225.0;
			double w2 = (322.0 + 13.0 * B) / 900.0;
			double w3 = (322.0 - 13.0 * B) / 900.0;
			double w_vec[5] = {w3, w2, w1, w2, w3};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 6)
		{
			double G_vec[6];
			double w_vec[6];

			G_vec[0] = -0.9324695142031520278123;
			G_vec[1] = -0.661209386466264513661;
			G_vec[2] = -0.2386191860831969086305;
			G_vec[3] = 0.238619186083196908631;
			G_vec[4] = 0.661209386466264513661;
			G_vec[5] = 0.9324695142031520278123;

			w_vec[0] = 0.1713244923791703450403;
			w_vec[1] = 0.3607615730481386075698;
			w_vec[2] = 0.4679139345726910473899;
			w_vec[3] = 0.46791393457269104739;
			w_vec[4] = 0.3607615730481386075698;
			w_vec[5] = 0.1713244923791703450403;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 8)
		{
			double G_vec[8];
			double w_vec[8];

			G_vec[0] = -0.9602898564975362316836;
			G_vec[1] = -0.7966664774136267395916;
			G_vec[2] = -0.5255324099163289858177;
			G_vec[3] = -0.1834346424956498049395;
			G_vec[4] = 0.1834346424956498049395;
			G_vec[5] = 0.5255324099163289858177;
			G_vec[6] = 0.7966664774136267395916;
			G_vec[7] = 0.9602898564975362316836;

			w_vec[0] = 0.1012285362903762591525;
			w_vec[1] = 0.2223810344533744705444;
			w_vec[2] = 0.313706645877887287338;
			w_vec[3] = 0.3626837833783619829652;
			w_vec[4] = 0.3626837833783619829652;
			w_vec[5] = 0.313706645877887287338;
			w_vec[6] = 0.222381034453374470544;
			w_vec[7] = 0.1012285362903762591525;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 10)
		{
			double G_vec[10];
			double w_vec[10];

			G_vec[0] = -0.9739065285171717;
			G_vec[1] = -0.8650633666889845;
			G_vec[2] = -0.6794095682990244;
			G_vec[3] = -0.4333953941292472;
			G_vec[4] = -0.1488743389816312;
			G_vec[5] = 0.1488743389816312;
			G_vec[6] = 0.4333953941292472;
			G_vec[7] = 0.6794095682990244;
			G_vec[8] = 0.8650633666889845;
			G_vec[9] = 0.9739065285171717;

			w_vec[0] = 0.0666713443086881;
			w_vec[1] = 0.1494513491505804;
			w_vec[2] = 0.2190863625159820;
			w_vec[3] = 0.2692667193099965;
			w_vec[4] = 0.2955242247147530;
			w_vec[5] = 0.2955242247147530;
			w_vec[6] = 0.2692667193099965;
			w_vec[7] = 0.2190863625159820;
			w_vec[8] = 0.1494513491505804;
			w_vec[9] = 0.0666713443086881;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 20)
		{
			double G_vec[20];
			double w_vec[20];

			G_vec[0]  = -0.99312859918509492478612239;
			G_vec[1]  = -0.96397192727791379126766613;
			G_vec[2]  = -0.91223442825132590586775244;
			G_vec[3]  = -0.83911697182221882339452906;
			G_vec[4]  = -0.74633190646015079261430507;
			G_vec[5]  = -0.6360536807265150254528367;
			G_vec[6]  = -0.51086700195082709800436405;
			G_vec[7]  = -0.3737060887154195606725482;
			G_vec[8]  = -0.2277858511416450780804962;
			G_vec[9]  = -0.07652652113349733375464041;
			G_vec[10] = 0.07652652113349733375464041;
			G_vec[11] = 0.2277858511416450780804962;
			G_vec[12] = 0.37370608871541956067254818;
			G_vec[13] = 0.5108670019508270980043641;
			G_vec[14] = 0.6360536807265150254528367;
			G_vec[15] = 0.74633190646015079261430507;
			G_vec[16] = 0.83911697182221882339452906;
			G_vec[17] = 0.91223442825132590586775244;
			G_vec[18] = 0.9639719272779137912676661;
			G_vec[19] = 0.9931285991850949247861224;

			w_vec[0]  = 0.01761400713915211831186196;
			w_vec[1]  = 0.04060142980038694133103995;
			w_vec[2]  = 0.0626720483341090635695065;
			w_vec[3]  = 0.083276741576704748724758143;
			w_vec[4]  = 0.10193011981724043503675014;
			w_vec[5]  = 0.1181945319615184173123774;
			w_vec[6]  = 0.1316886384491766268984945;
			w_vec[7]  = 0.14209610931838205132929833;
			w_vec[8]  = 0.14917298647260374678782874;
			w_vec[9]  = 0.15275338713072585069808433;
			w_vec[10] = 0.1527533871307258506980843;
			w_vec[11] = 0.1491729864726037467878287;
			w_vec[12] = 0.14209610931838205132929833;
			w_vec[13] = 0.1316886384491766268984945;
			w_vec[14] = 0.11819453196151841731237738;
			w_vec[15] = 0.1019301198172404350367501;
			w_vec[16] = 0.08327674157670474872475814;
			w_vec[17] = 0.0626720483341090635695065;
			w_vec[18] = 0.04060142980038694133103995;
			w_vec[19] = 0.01761400713915211831186196;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					w_ptr[temp_GP_1D * i + j] = w_vec[i] * w_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 0] = G_vec[j];
					point_ptr[(temp_GP_1D * i + j) * info->DIMENSION + 1] = G_vec[i];
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else
		{
			printf("Error, incorrect input data at gauss points.\nChange the number of the gauss points.\n");
			exit(1);
		}
	}
	else if (info->DIMENSION == 3)
	{
		int temp_GP_2D = temp_GP_1D * temp_GP_1D;
		if (temp_GP_1D == 3)
		{
			double G1 = pow((3.0 / 5.0), 0.5);
			double G_vec[3] = {-G1, 0.0, G1};
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;
			double w_vec[3] = {w2, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 4)
		{
			double A = pow((6.0 / 5.0), 0.5);
			double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
			double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
			double G_vec[4] = {-G2, -G1, G1, G2};
			double B = pow(30.0, 0.5);
			double w1 = (18.0 + B) / 36.0;
			double w2 = (18.0 - B) / 36.0;
			double w_vec[4] = {w2, w1, w1, w2};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 5)
		{
			double A = pow((10.0 / 7.0), 0.5);
			double G1 = pow((5.0 - 2.0 * A), 0.5) / 3.0;
			double G2 = pow((5.0 + 2.0 * A), 0.5) / 3.0;
			double G_vec[5] = {-G2, -G1, 0.0, G1, G2};
			double B = pow(70.0, 0.5);
			double w1 = 128.0 / 225.0;
			double w2 = (322.0 + 13.0 * B) / 900.0;
			double w3 = (322.0 - 13.0 * B) / 900.0;
			double w_vec[5] = {w3, w2, w1, w2, w3};

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 6)
		{
			double G_vec[6];
			double w_vec[6];

			G_vec[0] = -0.9324695142031520278123;
			G_vec[1] = -0.661209386466264513661;
			G_vec[2] = -0.2386191860831969086305;
			G_vec[3] = 0.238619186083196908631;
			G_vec[4] = 0.661209386466264513661;
			G_vec[5] = 0.9324695142031520278123;

			w_vec[0] = 0.1713244923791703450403;
			w_vec[1] = 0.3607615730481386075698;
			w_vec[2] = 0.4679139345726910473899;
			w_vec[3] = 0.46791393457269104739;
			w_vec[4] = 0.3607615730481386075698;
			w_vec[5] = 0.1713244923791703450403;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 8)
		{
			double G_vec[8];
			double w_vec[8];

			G_vec[0] = -0.9602898564975362316836;
			G_vec[1] = -0.7966664774136267395916;
			G_vec[2] = -0.5255324099163289858177;
			G_vec[3] = -0.1834346424956498049395;
			G_vec[4] = 0.1834346424956498049395;
			G_vec[5] = 0.5255324099163289858177;
			G_vec[6] = 0.7966664774136267395916;
			G_vec[7] = 0.9602898564975362316836;

			w_vec[0] = 0.1012285362903762591525;
			w_vec[1] = 0.2223810344533744705444;
			w_vec[2] = 0.313706645877887287338;
			w_vec[3] = 0.3626837833783619829652;
			w_vec[4] = 0.3626837833783619829652;
			w_vec[5] = 0.313706645877887287338;
			w_vec[6] = 0.222381034453374470544;
			w_vec[7] = 0.1012285362903762591525;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 10)
		{
			double G_vec[10];
			double w_vec[10];

			G_vec[0] = -0.9739065285171717;
			G_vec[1] = -0.8650633666889845;
			G_vec[2] = -0.6794095682990244;
			G_vec[3] = -0.4333953941292472;
			G_vec[4] = -0.1488743389816312;
			G_vec[5] = 0.1488743389816312;
			G_vec[6] = 0.4333953941292472;
			G_vec[7] = 0.6794095682990244;
			G_vec[8] = 0.8650633666889845;
			G_vec[9] = 0.9739065285171717;

			w_vec[0] = 0.0666713443086881;
			w_vec[1] = 0.1494513491505804;
			w_vec[2] = 0.2190863625159820;
			w_vec[3] = 0.2692667193099965;
			w_vec[4] = 0.2955242247147530;
			w_vec[5] = 0.2955242247147530;
			w_vec[6] = 0.2692667193099965;
			w_vec[7] = 0.2190863625159820;
			w_vec[8] = 0.1494513491505804;
			w_vec[9] = 0.0666713443086881;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else if (temp_GP_1D == 20)
		{
			double G_vec[20];
			double w_vec[20];

			G_vec[0]  = -0.99312859918509492478612239;
			G_vec[1]  = -0.96397192727791379126766613;
			G_vec[2]  = -0.91223442825132590586775244;
			G_vec[3]  = -0.83911697182221882339452906;
			G_vec[4]  = -0.74633190646015079261430507;
			G_vec[5]  = -0.6360536807265150254528367;
			G_vec[6]  = -0.51086700195082709800436405;
			G_vec[7]  = -0.3737060887154195606725482;
			G_vec[8]  = -0.2277858511416450780804962;
			G_vec[9]  = -0.07652652113349733375464041;
			G_vec[10] = 0.07652652113349733375464041;
			G_vec[11] = 0.2277858511416450780804962;
			G_vec[12] = 0.37370608871541956067254818;
			G_vec[13] = 0.5108670019508270980043641;
			G_vec[14] = 0.6360536807265150254528367;
			G_vec[15] = 0.74633190646015079261430507;
			G_vec[16] = 0.83911697182221882339452906;
			G_vec[17] = 0.91223442825132590586775244;
			G_vec[18] = 0.9639719272779137912676661;
			G_vec[19] = 0.9931285991850949247861224;

			w_vec[0]  = 0.01761400713915211831186196;
			w_vec[1]  = 0.04060142980038694133103995;
			w_vec[2]  = 0.0626720483341090635695065;
			w_vec[3]  = 0.083276741576704748724758143;
			w_vec[4]  = 0.10193011981724043503675014;
			w_vec[5]  = 0.1181945319615184173123774;
			w_vec[6]  = 0.1316886384491766268984945;
			w_vec[7]  = 0.14209610931838205132929833;
			w_vec[8]  = 0.14917298647260374678782874;
			w_vec[9]  = 0.15275338713072585069808433;
			w_vec[10] = 0.1527533871307258506980843;
			w_vec[11] = 0.1491729864726037467878287;
			w_vec[12] = 0.14209610931838205132929833;
			w_vec[13] = 0.1316886384491766268984945;
			w_vec[14] = 0.11819453196151841731237738;
			w_vec[15] = 0.1019301198172404350367501;
			w_vec[16] = 0.08327674157670474872475814;
			w_vec[17] = 0.0626720483341090635695065;
			w_vec[18] = 0.04060142980038694133103995;
			w_vec[19] = 0.01761400713915211831186196;

			for (int i = 0; i < temp_GP_1D; i++)
			{
				for (int j = 0; j < temp_GP_1D; j++)
				{
					for (int k = 0; k < temp_GP_1D; k++)
					{
						w_ptr[i * temp_GP_2D + j * temp_GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 0] = G_vec[k];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 1] = G_vec[j];
						point_ptr[(i * temp_GP_2D + j * temp_GP_1D + k) * info->DIMENSION + 2] = G_vec[i];
					}
				}
				w_1D_ptr[i] = w_vec[i];
				point_1D_ptr[i] = G_vec[i];
			}
		}
		else
		{
			printf("Error, incorrect input data at gauss points.\nChange the number of the gauss points.\n");
			exit(1);
		}
	}
}


void Make_Gauss_points(bool isSinglePatch, information *info)
{
	int n[2] = {0};
	n[0] = pow_int(info->c.NUM_GAUSS_POINTS, info->DIMENSION), n[1] = pow_int(info->c.NUM_GAUSS_POINTS_EXTENDED, info->DIMENSION);

	Make_gauss_array(info);
	if (isSinglePatch)
	{
		bool opp_flag = false;
		#pragma omp parallel for
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];
			info->gp[i].setVar(opp_flag, p, i, n[0], info);
		}
	}
	else
	{
		bool opp_flag = true;
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];

			// normal gauss points
			if (info->c.USE_EXTENDED_QUADRATURE == 0)
				info->gp[i].setVar(opp_flag, p, i, n[0], info);

			// extended gauss points
			else if (info->c.USE_EXTENDED_QUADRATURE == 1)
			{
				if (info->eoi[i].size() > 1)
				{
					if (i < info->Total_Element_to_mesh[1])
						info->gp[i].setVar(opp_flag, p, i, n[0], info);
					else
						info->gp[i].setVar(opp_flag, p, i, n[0], n[1], info);
				}
				else
					info->gp[i].setVar(opp_flag, p, i, n[0], info);
			}

			// error message
			else
			{
				printf("info->c.USE_EXTENDED_QUADRATURE is not correct.\n");
				exit(1);
			}
		}
	}
}


void gp_switch(bool flag, information *info)
{
	// if flag == true, switch to the extended gauss points. if flag == false, switch to the original gauss points.
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		info->gp[i].switch_gp(flag);
	}
}


// First Inverese mapping: ニュートン法によってグローバルIGA解析モデルの物理座標から自然座標上のモデルを探索する
void First_inverse_mapping(information *info)
{
	int i,j,k;
	double MAX_ITER = 1000;
	double MAX_RESIDUAL = 1.0e-10;
	double residual_norm = 0.0;
	vector<double> para_patch_coord(info->DIMENSION, 0.5); // 初期値は自然座標の中心点

	int vertex_num;
	if (info->DIMENSION == 2)
		vertex_num = 9;
	else if (info->DIMENSION == 3)
		vertex_num = 27;

	int local_ele_num = info->Geo_Total_Element_on_mesh;
	
	for (i = 0; i < local_ele_num; i++)
	{
		for (j = 0; j < vertex_num; j++)
		{
			for (k = 0; k < MAX_ITER; k++)
			{
				// calculate residulal vector: R(ξ) = N(ξ) * X - x
				vector<double> Residual_vector(info->DIMENSION);
				vector<double> para_ele_coord(info->DIMENSION, 0.0);
				vector<double> pc(info->DIMENSION, 0.0);

				int temp_ele = ele_check(info->Global_local_patch, para_patch_coord.data(), info);
				tilde_coord(para_ele_coord.data(), para_patch_coord.data(), info->Element_patch[temp_ele], temp_ele,info);
				physical_coord(temp_ele, para_ele_coord.data(), pc.data(), info);
				for (int d = 0; d < info->DIMENSION; d++)
					Residual_vector[d] = pc[d] - info->Target_Physical_Coord[(i * vertex_num + j) * info->DIMENSION + d];

				// compute norm (reset each iteration)
				residual_norm = 0.0;
				for (int d = 0; d < info->DIMENSION; d++)
					residual_norm += Residual_vector[d] * Residual_vector[d];
				residual_norm = sqrt(residual_norm);

				// debug
				printf("Element %d, Vertex %d, Iteration %d: para_patch_coord = (", i, j, k);
				for (int d = 0; d < info->DIMENSION; d++)
				    printf("%f ", para_patch_coord[d]);
				printf(")\n");

				// printf("Element %d, Vertex %d, Iteration %d: Residual Norm = %e\n", i, j, k, residual_norm);

				// convergence check
				if (residual_norm < MAX_RESIDUAL)
				{
					for (int d = 0; d < info->DIMENSION; d++)
						info->Target_Para_Coord[(i * vertex_num + j) * info->DIMENSION + d] = para_patch_coord[d];
					break;
				}
				
				// calculate delta ξ: -(J)^{-1} * R(ξ)
				vector<double> delta_para(info->DIMENSION, 0.0);
				vector<double> Jacobian_inverse(info->DIMENSION * info->DIMENSION);
				vector<double> dummy_B(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
				Make_B_component(temp_ele, para_ele_coord.data(), dummy_B.data(), info, 1, Jacobian_inverse.data());

				// Note: Jacobian is stored in column-major order (index = col*DIM + row)
				for (int d = 0; d < info->DIMENSION; d++) {
					double sum = 0.0;
					for (int e = 0; e < info->DIMENSION; e++)
						sum += Jacobian_inverse[e * info->DIMENSION + d] * Residual_vector[e];
					delta_para[d] = -sum;
				}
				
				// update ξ: ξ = ξ + delta ξ
				for (int d = 0; d < info->DIMENSION; d++)
					para_patch_coord[d] += delta_para[d];

				// ξは0.0~1.0の間しかとらない。超える場合は最小、または最大値に修正
				for (int d = 0; d < info->DIMENSION; d++)
				{
					if (para_patch_coord[d] >= 1.0)
						para_patch_coord[d] = 1.0;
					if (para_patch_coord[d] <= 0.0)
						para_patch_coord[d] = 0.0;
				}
			}
			if (k == MAX_ITER)
			{
				printf("ERROR: Maximum iterations reached without convergence for vertex %d of element %d\n", j, i);
				exit(1);
			}
			else 
				printf("Vertex %d of element %d converged in %d iterations with residual norm %e\n", j, i, k, residual_norm);
		}
	}
	// debug
	for (i = 0; i < local_ele_num; i++)
	{
		for (j = 0; j < vertex_num; j++)
		{
			printf("Element %d, Vertex %d: Target Physical Coord = (", i, j);
			for (int d = 0; d < info->DIMENSION; d++)
				printf("%f ", info->Target_Physical_Coord[(i * vertex_num + j) * info->DIMENSION + d]);
			printf("), Target Para Coord = (");
			for (int d = 0; d < info->DIMENSION; d++)
				printf("%f ", info->Target_Para_Coord[(i * vertex_num + j) * info->DIMENSION + d]);
			printf(")\n");
		}
	}
}


// Second Inverse mapping: 最小二乗法によって、ローカルIGA解析モデルの制御点を移動させる
void Second_inverse_mapping(information *info)
{
	const int fasten_flag_x = 0; // 0: fasten, 1: calculate
	const int fasten_flag_y = 0; // 0: fasten, 1: calculate
	const int fasten_flag_z = 1; // 0: fasten, 1: calculate
	const int element_num = info->Geo_Total_Element_on_mesh;
	const int sample_num = (info->DIMENSION == 2) ? 9 : 27;

	for (int cp = 0; cp < info->Geo_Total_Control_Point_on_mesh; cp++)
	{
		for (int d = 0; d < info->DIMENSION + 1; d++)
			info->New_Node_Coordinate[cp * (info->DIMENSION + 1) + d] = info->Geo_Node_Coordinate[cp * (info->DIMENSION + 1) + d];
	}

	for (int ele = 0; ele < element_num; ele++)
	{
		const int patch = info->Element_patch[ele];
		const int cp_num = info->No_Control_point_ON_ELEMENT[patch];

		vector<double> shape_matrix(sample_num * cp_num, 0.0);
		vector<double> target(sample_num * info->DIMENSION, 0.0);

		for (int j = 0; j < sample_num; j++)
		{
			vector<double> vertex_para(info->DIMENSION, 0.0);
			if (info->DIMENSION == 2)
			{
				const double sample_coord[9][2] = {
					{-1.0, -1.0}, {0.0, -1.0}, {1.0, -1.0},
					{-1.0,  0.0}, {0.0,  0.0}, {1.0,  0.0},
					{-1.0,  1.0}, {0.0,  1.0}, {1.0,  1.0}
				};
				vertex_para[0] = sample_coord[j][0];
				vertex_para[1] = sample_coord[j][1];
			}
			else if (info->DIMENSION == 3)
			{
				const double sample_coord[27][3] = {
					{-1.0, -1.0, -1.0}, { 0.0, -1.0, -1.0}, { 1.0, -1.0, -1.0},
					{-1.0,  0.0, -1.0}, { 0.0,  0.0, -1.0}, { 1.0,  0.0, -1.0},
					{-1.0,  1.0, -1.0}, { 0.0,  1.0, -1.0}, { 1.0,  1.0, -1.0},

					{-1.0, -1.0,  0.0}, { 0.0, -1.0,  0.0}, { 1.0, -1.0,  0.0},
					{-1.0,  0.0,  0.0}, { 0.0,  0.0,  0.0}, { 1.0,  0.0,  0.0},
					{-1.0,  1.0,  0.0}, { 0.0,  1.0,  0.0}, { 1.0,  1.0,  0.0},

					{-1.0, -1.0,  1.0}, { 0.0, -1.0,  1.0}, { 1.0, -1.0,  1.0},
					{-1.0,  0.0,  1.0}, { 0.0,  0.0,  1.0}, { 1.0,  0.0,  1.0},
					{-1.0,  1.0,  1.0}, { 0.0,  1.0,  1.0}, { 1.0,  1.0,  1.0}
				};
				vertex_para[0] = sample_coord[j][0];
				vertex_para[1] = sample_coord[j][1];
				vertex_para[2] = sample_coord[j][2];
			}

			vector<double> R(MAX_NO_CP_ON_ELEMENT, 0.0);
			vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION, 0.0);
			geo_shape_and_dshape(R.data(), dR.data(), vertex_para.data(), ele, true, info);

			for (int cp = 0; cp < cp_num; cp++)
				shape_matrix[j * cp_num + cp] = R[cp];

			for (int d = 0; d < info->DIMENSION; d++)
				target[j * info->DIMENSION + d] = info->Target_Para_Coord[(ele * sample_num + j) * info->DIMENSION + d];
		}

		// Solve normal equations: (A^T * A) * x = A^T * b
		for (int d = 0; d < info->DIMENSION; d++)
		{
			if ((d == 0 && fasten_flag_x == 0) || (d == 1 && fasten_flag_y == 0) || (d == 2 && fasten_flag_z == 0))
				continue;

			vector<double> normal_matrix(cp_num * cp_num, 0.0);
			vector<double> rhs(cp_num, 0.0);
			vector<double> solution(cp_num, 0.0);

			for (int cp1 = 0; cp1 < cp_num; cp1++)
			{
				for (int cp2 = 0; cp2 < cp_num; cp2++)
				{
					double sum = 0.0;
					for (int v = 0; v < sample_num; v++)
						sum += shape_matrix[v * cp_num + cp1] * shape_matrix[v * cp_num + cp2];
					normal_matrix[cp1 * cp_num + cp2] = sum;
				}

				double sum = 0.0;
					for (int v = 0; v < sample_num; v++)
					sum += shape_matrix[v * cp_num + cp1] * target[v * info->DIMENSION + d];
				rhs[cp1] = sum;
			}

			GaussianElimination(solution.data(), rhs.data(), normal_matrix.data(), cp_num);

			for (int cp = 0; cp < cp_num; cp++)
			{
				const int global_cp = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + cp];
				info->New_Node_Coordinate[global_cp * (info->DIMENSION + 1) + d] = solution[cp];
			}
		}
	}
}

// calcurate B
double Make_B_component(int ele, double *para, double *out_B_component, information *info, int return_mode, double *out_a_matrix)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];

	vector<double> R(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	shape_and_dshape(R.data(), dR.data(), para, ele, true, info);

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					a_2x2[i * 2 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id];
				}
			}
		}

		if (return_mode == 2 && out_a_matrix != nullptr)
		{
			for (int i = 0; i < 4; i++)
				out_a_matrix[i] = a_2x2[i];
		}

		J = InverseMatrix_2x2(a_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k * 2 + i] * dR[j * info->DIMENSION + k];
			}
		}

		// Copy inverse a_2x2 to out_a_matrix if return_mode == 1
		if (return_mode == 1 && out_a_matrix != nullptr)
		{
			for (int i = 0; i < 4; i++)
				out_a_matrix[i] = a_2x2[i];
		}

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; k++)
				{
					int id = info->Controlpoint_of_Element[ele * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					a_3x3[i * 3 + j] += dR[k * info->DIMENSION + j] * info->Node_Coordinate[id];
				}
			}
		}

		if (return_mode == 2 && out_a_matrix != nullptr)
		{
			for (int i = 0; i < 9; i++)
				out_a_matrix[i] = a_3x3[i];
		}

		J = InverseMatrix_3x3(a_3x3);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k * 3 + i] * dR[j * info->DIMENSION + k];
			}
		}

		// Copy inverse a_3x3 to out_a_matrix if return_mode == 1
		if (return_mode == 1 && out_a_matrix != nullptr)
		{
			for (int i = 0; i < 9; i++)
				out_a_matrix[i] = a_3x3[i];
		}

		return J;
	}

	return J;
}



// SSIGA用のMake_B_component
double Make_B_component_for_SSIGA(int ele_disp, double *para_disp, double *out_B_component, information *info, int return_mode, double *out_c_matrix)
{
	double J = 0.0;
	double a_2x2[4], a_3x3[9];
	double b_2x2[4], b_3x3[9];
	double c_2x2[4], c_3x3[9];

	// パラメータ座標再定義　(para_dispはローカル変位要素のパラメータ座標)
	vector<double> para_geo(info->DIMENSION); // ローカル形状要素パラメータ座標
	vector<double> para_glo(info->DIMENSION); // グローバル要素パラメータ座標

	// 要素番号の取得と座標系の変換
	int ele_geo = trans_local_para_to_local_geo_para(ele_disp, para_disp, para_geo.data(), info);
	int ele_glo = trans_local_para_to_global_para(ele_disp, para_disp, para_glo.data(), info);

	vector<double> R_disp(MAX_NO_CP_ON_ELEMENT);
	vector<double> R_geo(MAX_NO_CP_ON_ELEMENT);
	vector<double> R_glo(MAX_NO_CP_ON_ELEMENT);
	vector<double> dR_disp(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	vector<double> dR_geo(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	vector<double> dR_glo(MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	Bspline_shape_and_dshape(R_disp.data(), dR_disp.data(), para_disp, ele_disp, true, info);
	geo_shape_and_dshape(R_geo.data(), dR_geo.data(), para_geo.data(), ele_geo, false, info);
	shape_and_dshape(R_glo.data(), dR_glo.data(), para_glo.data(), ele_glo, false, info);

	
	// 変位要素空間の微分によって生じる定数
	for (int i = 0; i < info->DIMENSION; i++)
	{
		double coef = dShapeFunc_from_paren(i, ele_disp, info);
		for (int j = 0; j < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; j++)
			dR_geo[j * info->DIMENSION + i] *= coef;
	}

	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; k++)
				{
					int id = info->Geo_Controlpoint_of_Element[ele_geo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					a_2x2[i * 2 + j] += dR_geo[k * info->DIMENSION + j] * info->Geo_Node_Coordinate[id];
				}
			}
		}
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				b_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_glo]]; k++)
				{
					int id = info->Controlpoint_of_Element[ele_glo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					b_2x2[i * 2 + j] += dR_glo[k * info->DIMENSION + j] * info->Node_Coordinate[id];
				}
			}
		}
        for (int i = 0; i < 2; i++)
        {
			for (int j = 0; j < 2; j++)
            {
				c_2x2[i * 2 + j] = 0.0;
				for (int k = 0; k < 2; k++)
                {
					c_2x2[i * 2 + j] += b_2x2[i * 2 + k] * a_2x2[k * 2 + j];
                }
            }
        }

		if (return_mode == 2 && out_c_matrix != nullptr)
		{
			for (int i = 0; i < 4; i++)
				out_c_matrix[i] = c_2x2[i];
		}

		J = InverseMatrix_2x2(c_2x2);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_disp]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += c_2x2[k * 2 + i] * dR_disp[j * info->DIMENSION + k];
			}
		}

		// Copy inverse c_2x2 to out_c_matrix if return_mode == 1
		if (return_mode == 1 && out_c_matrix != nullptr)
		{
			for (int i = 0; i < 4; i++)
				out_c_matrix[i] = c_2x2[i];
		}

		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->Geo_No_Control_point_ON_ELEMENT[info->Geo_Element_patch[ele_geo]]; k++)
				{
					int id = info->Geo_Controlpoint_of_Element[ele_geo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					a_3x3[i * 3 + j] += dR_geo[k * info->DIMENSION + j] * info->Geo_Node_Coordinate[id];
				}
			}
		}
		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->DIMENSION; j++)
			{
				b_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_glo]]; k++)
				{
					int id = info->Controlpoint_of_Element[ele_glo * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i;
					b_3x3[i * 3 + j] += dR_glo[k * info->DIMENSION + j] * info->Node_Coordinate[id];
				}
			}
		}
        for (int i = 0; i < 3; i++)
        {
			for (int j = 0; j < 3; j++)
            {
				c_3x3[i * 3 + j] = 0.0;
				for (int k = 0; k < 3; k++)
                {
                    c_3x3[i * 3 + j] += b_3x3[i * 3 + k] * a_3x3[k * 3 + j];
                }
            }
        }

		if (return_mode == 2 && out_c_matrix != nullptr)
		{
			for (int i = 0; i < 9; i++)
				out_c_matrix[i] = c_3x3[i];
		}

		J = InverseMatrix_3x3(c_3x3);

		for (int i = 0; i < info->DIMENSION; i++)
		{
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele_disp]]; j++)
			{
				out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (int k = 0; k < info->DIMENSION; k++)
					out_B_component[i * MAX_NO_CP_ON_ELEMENT + j] += c_3x3[k * 3 + i] * dR_disp[j * info->DIMENSION + k];
			}
		}

		// Copy inverse c_3x3 to out_c_matrix if return_mode == 1
		if (return_mode == 1 && out_c_matrix != nullptr)
		{
			for (int i = 0; i < 9; i++)
				out_c_matrix[i] = c_3x3[i];
		}

		return J;
	}

	return J;
}


int trans_local_para_to_local_geo_para(int ele_disp, double *para_disp, double *para_geo, information *info)
{
	int i;
	int ele_geo;

	int disp_patch_num = info->Element_patch[ele_disp];		//要素番号からパッチ番号を取得
	int geo_patch_num = disp_patch_num - info->Total_Patch_on_mesh[0];	//ジオメトリ表現のパッチ番号は、ローカルのパッチ番号と常に等しい

	double para_coord_tilde[MAX_DIMENSION] = {0.0};			//ローカルの要素パラメータ座標
	double para_coord[MAX_DIMENSION] = {0.0};				//ローカルのパッチパラメータ座標
	double para_coord_tilde_geo[MAX_DIMENSION] = {0.0};		//ローカルジオメトリ表現の要素パラメータ座標

	for (i = 0; i < info->DIMENSION; i++)
		para_coord_tilde[i] = para_disp[i];							//要素のパラメータ座標をガウス点座標に変換
	trans_ele_patch_coord(para_coord, para_coord_tilde, disp_patch_num, ele_disp, info);	//要素のパラメータ座標をパッチのパラメータ座標に変換
	ele_geo = geo_ele_check(geo_patch_num, para_coord, info);							//パラメータ座標から、要素を探索
	geo_tilde_coord(para_coord_tilde_geo, para_coord, geo_patch_num, ele_geo, info);			//パッチパラメータ座標を、ジオメトリ要素パラメータ座標に変換

	for (i = 0; i < info->DIMENSION; i++)
		para_geo[i] = para_coord_tilde_geo[i];	//要素のパラメータ座標をガウス点座標に変換

	return ele_geo;
}


int trans_local_para_to_global_para(int ele_loc, double *para_loc, double *para_glo, information *info)
{
	int i;
	int ele_glo;
	int local_patch_num = info->Element_patch[ele_loc];		//要素番号からパッチ番号を取得
	int geo_patch_num = local_patch_num - info->Total_Patch_on_mesh[0];	//ジオメトリ表現のパッチ番号は、ローカルのパッチ番号と常に等しい

	double para_coord_glo[MAX_DIMENSION] = {0.0};			//グローバルのパラメータ空間座標
	double para_coord_tilde[MAX_DIMENSION] = {0.0};			//ローカルの要素パラメータ座標
	double para_coord[MAX_DIMENSION] = {0.0};				//ローカルのパッチパラメータ座標
	double para_geo_coord_tilde[MAX_DIMENSION] = {0.0};		//ローカルジオメトリ表現の要素パラメータ座標
	for (i = 0; i < info->DIMENSION; i++)
		para_coord_tilde[i] = para_loc[i];													//要素のパラメータ座標をガウス点座標に変換
	trans_ele_patch_coord(para_coord, para_coord_tilde, local_patch_num, ele_loc, info);	//要素のパラメータ座標をパッチのパラメータ座標に変換
	int geo_e = geo_ele_check(geo_patch_num, para_coord, info);								//パラメータ座標から、要素を探索
	geo_tilde_coord(para_geo_coord_tilde, para_coord, geo_patch_num, geo_e, info);			//パッチパラメータ座標を、ジオメトリ要素パラメータ座標に変換

	geo_parameter_coord(geo_e, para_geo_coord_tilde, para_coord_glo, info);					//ジオメトリ要素パラメータ座標を、グローバルのパラメータ空間座標に変換

	ele_glo = ele_check(info->Global_local_patch, para_coord_glo, info); 					// グローバルのパラメータ空間座標から要素番号を計算	
	tilde_coord(para_glo, para_coord_glo, info->Global_local_patch, ele_glo, info);
	
	return ele_glo;
}


double Make_updated_jacobian(int ele, double *para, information *info)
{
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	double J = Make_B_component(ele, para, b.data(), info);
	return J;
}


double Make_B_Linear(int ele, double *para, double *out_BL, information *info)
{
	double J = 0.0;
	vector<double> b(info->DIMENSION * MAX_NO_CP_ON_ELEMENT);
	if (ele < info->Total_Element_to_mesh[1])
		J = Make_B_component(ele, para, b.data(), info);
	else
		J = Make_B_component_for_SSIGA(ele, para, b.data(), info);
	if (info->DIMENSION == 2)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	else if (info->DIMENSION == 3)
	{
		for (int i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[ele]]; i++)
		{
			// Linear
			out_BL[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			out_BL[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
		return J;
	}
	return J;
}


// K matrix
void Make_D_Matrix(information *info)
{
	int i, j;

	if (info->DIMENSION == 2)
	{
		if (info->c.ANALYSIS_MODE == 0) // 平面応力状態
		{
			double Eone = E / (1.0 - nu * nu);
			double D1[3][3] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) * Eone / 2.0}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
		else if (info->c.ANALYSIS_MODE == 1) // 平面ひずみ状態
		{
			double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			double D1[3][3] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2.0 / (1.0 - nu) * Eone}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
	}
	else if (info->DIMENSION == 3)
	{
		double E_ii = (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E;
		double E_ij = nu / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E;
		double G = E / (2.0 * (1.0 + nu));

		double D1[6][6] = {{E_ii, E_ij, E_ij, 0.0, 0.0, 0.0},
						   {E_ij, E_ii, E_ij, 0.0, 0.0, 0.0},
						   {E_ij, E_ij, E_ii, 0.0, 0.0, 0.0},
						   { 0.0, 0.0,  0.0,  G  , 0.0, 0.0},
						   { 0.0, 0.0,  0.0,  0.0, G  , 0.0},
						   { 0.0, 0.0,  0.0,  0.0, 0.0, G  }};

		for (i = 0; i < D_MATRIX_SIZE; i++)
			for (j = 0; j < D_MATRIX_SIZE; j++)
				info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
	}
}


// 拘束されている行数を省いた行列の番号の制作
void Make_Index_Dof(information *info)
{
	int i, k = 0;

	// 拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		info->Index_Dof[info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1]] = ERROR;
	}
	// ERROR以外に番号を付ける
	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		if (info->Index_Dof[i] != ERROR)
		{
			info->Index_Dof[i] = k;
			k++;
		}
	}

	K_Whole_Size = k;
	printf("\nK_Whole_Size = %d\n\n", k);
}


void Make_K_Whole_Ptr_Col(information *info, int mode_select)
{
	// mode_select == 0: Ptr を作成
	// mode_select == 1: Col を作成
	
	vector<set<int>> &ptr = info->ptr;

	if (mode_select == 0)
	{
		ptr.resize(K_Whole_Size);
		vector<std::mutex> mtx(K_Whole_Size);

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			int p = info->Element_patch[i];
			for (int j = 0; j < info->No_Control_point_ON_ELEMENT[p]; j++)
				for (int k = 0; k < info->DIMENSION; k++)
				{
					int row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					if (row >= 0)
						for (int l = 0; l < info->No_Control_point_ON_ELEMENT[p]; l++)
							for (int m = 0; m < info->DIMENSION; m++)
							{
								int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + m];
								if (col >= row)
								{
									std::lock_guard<std::mutex> lock(mtx[row]);
									ptr[row].insert(col);
								}
							}
				}

			// coupling element
			if (i >= info->Total_Element_to_mesh[1])
				for (size_t j = 0; j < info->eoi[i].size(); j++)
					for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; k++)
						for (int l = 0; l < info->DIMENSION; l++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + l];
							if (row >= 0)
								for (int m = 0; m < info->No_Control_point_ON_ELEMENT[p]; m++)
									for (int n = 0; n < info->DIMENSION; n++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + n];
										if (col >= row)
										{
											std::lock_guard<std::mutex> lock(mtx[row]);
											ptr[row].insert(col);
										}
									}
						}
		}

		// substitute ptr to K_Whole_Ptr
		info->K_Whole_Ptr[0] = 0;
		for (int i = 0; i < K_Whole_Size; i++)
			info->K_Whole_Ptr[i + 1] = info->K_Whole_Ptr[i] + ptr[i].size();
	}
	
	else if (mode_select == 1)
	{
		#pragma omp parallel for
		for (int i = 0; i < K_Whole_Size; i++)
		{
			int count = 0;
			for (auto itr = ptr[i].begin(); itr != ptr[i].end(); ++itr)
			{
				info->K_Whole_Col[info->K_Whole_Ptr[i] + count] = *itr;
				count++;
			}
		}
	}
}


// Print progress bar (using \r to overwrite)
void print_progress_bar(const char* label, int current, int total)
{
	const int BAR_WIDTH = 30;
	int percent = (100 * current) / total;
	int filled = (percent * BAR_WIDTH) / 100;
	
	fprintf(stdout, "\r[%s]: %3d%% ¦¦ ", label, percent);
	for (int i = 0; i < BAR_WIDTH; i++) {
		if (i < filled) fprintf(stdout, "█");
		else fprintf(stdout, " ");
	}
	fprintf(stdout, "¦¦ (%d/%d elements) ", current, total);
	fflush(stdout);
}

// K matrix の値を求める
void Calc_K(information *info)
{
	int total_elements = info->Total_Element_to_mesh[Total_mesh];
	int count_processed = 0;
	
	fprintf(stdout, "Calculating K matrix...\n");
	fflush(stdout);
	
	gp_switch(false, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		vector<double> K_linear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		Calc_K_linear_EL(i, K_linear.data(), info);

		// K_EL add to K_Whole_Val
		for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
		{
			for (int j2 = 0; j2 < info->DIMENSION; j2++)
			{
				long row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
				if (row >= 0)
					for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
						for (int k2 = 0; k2 < info->DIMENSION; k2++)
						{
							long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
							if (col >= row)
							{
								int id = RowCol_to_icount(row, col, info);
								double val = K_linear[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];
								if (id < 0) continue;

								#pragma omp atomic
								info->K_Whole_Val[id] += val;
							}
						}
			}
		}

		
		// Progress bar output (every 1% or when completed)
		#pragma omp critical
		{
			count_processed++;
			int progress_interval = max(1, total_elements / 100);
			if (count_processed % progress_interval == 0 || count_processed == total_elements)
			{
				print_progress_bar("K matrix", count_processed, total_elements);
			}
		}
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "K matrix calculation completed.\n");
	fflush(stdout);

	if (Total_mesh < 2)
		return;

	// for coupling
	fprintf(stdout, "Calculating coupled K matrix...\n");
	fflush(stdout);
	
	count_processed = 0;
	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < total_elements; i++)
	{
		vector<double> K_linear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		// coupling
		if (Total_mesh > 1 && (info->Element_mesh[i] > 0 && info->eoi[i].size() > 0))
		{
			for (size_t j = 0; j < info->eoi[i].size(); j++)
			{
				Calc_coupled_K_linear_EL(i, info->eoi[i][j], K_linear_coupled.data(), info);

				// coupled_K_EL add to K_Whole_Val
				for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
					for (int j2 = 0; j2 < info->DIMENSION; j2++)
					{
						long row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
						if (row >= 0)
							for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
								for (int k2 = 0; k2 < info->DIMENSION; k2++)
								{
									long col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
									if (col >= row)
									{
										int id = RowCol_to_icount(row, col, info);
										if (id < 0) continue;  // Guard against structure mismatch
										
										double val = K_linear_coupled[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];
										
										#pragma omp atomic
										info->K_Whole_Val[id] += val;
									}
								}
					}
			}
		}
		
		// Progress bar output for coupling
		#pragma omp critical
		{
			count_processed++;
			int progress_interval = max(1, total_elements / 100);
			if (count_processed % progress_interval == 0 || count_processed == total_elements)
			{
				print_progress_bar("Coupled K", count_processed, total_elements);
			}
		}
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "K matrix calculation completed.\n");
	fflush(stdout);
}


void Calc_K_linear_EL(int El_No, double *K_EL, information *info)
{
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]] * info->DIMENSION;
	vector<double> B_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No].n(); i++)
	{
		// Generate linear B matrix
		double J = Make_B_Linear(El_No, &info->gp[El_No].para()[i * info->DIMENSION], B_linear.data(), info);

		// calc K_EL
		BDBJ(KIEL_SIZE, B_linear.data(), info->D, J, K1.data());
		for (int j = 0; j < KIEL_SIZE; j++)
			for (int k = 0; k < KIEL_SIZE; k++)
				K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No].w()[i] * K1[j * MAX_KIEL_SIZE + k];
	}
}


void Calc_coupled_K_linear_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;
	int KIEL_SIZE_loc = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_loc]] * info->DIMENSION;

	vector<double> BL_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> BG_linear(D_MATRIX_SIZE * MAX_KIEL_SIZE);
	vector<double> K1(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;

	for (int i = 0; i < info->gp[El_No_loc].n(); i++)
	{
		if (!(info->gp[El_No_loc].isOverlay()[i] && (info->gp[El_No_loc].opp_ele()[i] == El_No_glo)))
			continue;

		// Generate local linear B matrix
		double J = Make_B_Linear(El_No_loc, &info->gp[El_No_loc].para()[i * info->DIMENSION], BL_linear.data(), info);

		// // Debugging
		// cout << "El_No_loc: " << El_No_loc << ", El_No_glo: " << El_No_glo << ", GP index: " << i << ", Para coord: ";
		// cout << info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION + 0] << ", "
		// 	 << info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION + 1] << endl;
		
		// Generate global linear B matrix
		Make_B_Linear(El_No_glo, &info->gp[El_No_loc].opp_para_tilde()[i * info->DIMENSION], BG_linear.data(), info);

		coupled_BDBJ(KIEL_SIZE_glo, KIEL_SIZE_loc, BL_linear.data(), BG_linear.data(), info->D, J, K1.data());
		for (int j = 0; j < KIEL_SIZE_glo; j++)
			for (int k = 0; k < KIEL_SIZE_loc; k++)
				coupled_K_EL[j * MAX_KIEL_SIZE + k] += info->gp[El_No_loc].w()[i] * K1[j * MAX_KIEL_SIZE + k];
	}
}


void BDBJ(int KIEL_SIZE, double *B, double *D, double J, double *K_EL)
{
	vector<double> BD(MAX_KIEL_SIZE * D_MATRIX_SIZE);

	for (int i = 0; i < KIEL_SIZE * D_MATRIX_SIZE; i++)
		BD[i] = 0.0;

	// [B]T[D][B]Jの計算
	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				BD[i * D_MATRIX_SIZE + j] += B[k * MAX_KIEL_SIZE + i] * D[k * D_MATRIX_SIZE + j];

	for (int i = 0; i < KIEL_SIZE; i++)
		for (int j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
}


void coupled_BDBJ(int KIEL_SIZE_glo, int KIEL_SIZE_loc, double *B, double *BG, double *D, double J, double *K_EL)
{
	vector<double> BD(MAX_KIEL_SIZE * D_MATRIX_SIZE);

	for (int i = 0; i < KIEL_SIZE_glo * D_MATRIX_SIZE; i++)
		BD[i] = 0.0;

	//[B]GT[D][B]LJの計算
	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < D_MATRIX_SIZE; j++)
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				BD[i * D_MATRIX_SIZE + j] += BG[k * MAX_KIEL_SIZE + i] * D[k * D_MATRIX_SIZE + j];

	for (int i = 0; i < KIEL_SIZE_glo; i++)
		for (int j = 0; j < KIEL_SIZE_loc; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (int k = 0; k < D_MATRIX_SIZE; k++)
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
}


// F vector
void Make_F_Vec(information *info)
{
	int i, index;

	for (i = 0; i < info->Total_Load_to_mesh[Total_mesh]; i++)
	{
		index = info->Index_Dof[info->Load_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Load_Node_Dir[i * 2 + 1]];
		if (index >= 0)
			info->rhs_vec[index] += info->Value_of_Load[i];
	}
}


// 強制変位
void Make_F_Vec_disp_const(information *info)
{
	gp_switch(false, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		vector<double> K_linear(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		int disp_const_counter = 0;
		for (int j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
			for (int k = 0; k < info->DIMENSION; k++)
			{
				int id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
				if (id < 0)
				{
					#pragma omp atomic
					disp_const_counter++;
				}
			}

		if (disp_const_counter > 0)
		{
			Calc_K_linear_EL(i, K_linear.data(), info);

			for (int j = 0; j < info->DIMENSION; j++)
			{
				for (int k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k++)
				{
					int ii = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k] * info->DIMENSION + j;
					int b = info->Index_Dof[ii];
					if (b >= 0)
					{
						int ii_local = k * info->DIMENSION + j;
						for (int l = 0; l < info->DIMENSION; l++)
						{
							for (int m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
							{
								int jj = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l;
								int bb = info->Index_Dof[jj];
								if (bb < 0)
								{
									int jj_local = m * info->DIMENSION + l;
									for (int kk_const = 0; kk_const < info->Total_Constraint_to_mesh[Total_mesh]; kk_const++)
									{
										if (info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] == info->Constraint_Node_Dir[kk_const * 2 + 0] && l == info->Constraint_Node_Dir[kk_const * 2 + 1])
										{
											double val = -K_linear[ii_local * MAX_KIEL_SIZE + jj_local] * info->Value_of_Constraint[kk_const];

											#pragma omp atomic
											info->rhs_vec[b] += val;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	gp_switch(true, info);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		vector<double> K_linear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);
		vector<double> K_nonlinear_coupled(MAX_KIEL_SIZE * MAX_KIEL_SIZE);

		if (Total_mesh >= 2 && (info->Element_mesh[i] > 0 && info->eoi[i].size() > 0))
		{
			for (size_t j = 0; j < info->eoi[i].size(); j++)
			{
				int loc_disp_const_counter = 0;
				for (int k = 0; k < info->DIMENSION; k++)
				{
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; l++)
					{
						int loc_id = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
						if (loc_id < 0)
						{
							#pragma omp atomic
							loc_disp_const_counter++;
						}
					}
				}

				int glo_disp_const_counter = 0;
				for (int k = 0; k < info->DIMENSION; k++)
				{
					for (int l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; l++)
					{
						int glo_id = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + l] * info->DIMENSION + k];
						if (glo_id < 0)
						{
							#pragma omp atomic
							glo_disp_const_counter++;
						}
					}
				}

				if (loc_disp_const_counter > 0 || glo_disp_const_counter > 0)
					Calc_coupled_K_linear_EL(i, info->eoi[i][j], K_linear_coupled.data(), info);

				if (loc_disp_const_counter > 0)
				{
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
							{
								int row_local = j1 * info->DIMENSION + j2;
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (col < 0)
										{
											int col_local = k1 * info->DIMENSION + k2;
											for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
											{
												if (info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
												{
													double val = -K_linear_coupled[row_local * MAX_KIEL_SIZE + col_local] * info->Value_of_Constraint[n];

													#pragma omp atomic
													info->rhs_vec[row] += val;
												}
											}
										}
									}
							}
						}
				}

				if (glo_disp_const_counter > 0)
				{
					for (int j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; j1++)
						for (int j2 = 0; j2 < info->DIMENSION; j2++)
						{
							int row = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
							if (row >= 0)
							{
								int row_local = j1 * info->DIMENSION + j2;
								for (int k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->eoi[i][j]]]; k1++)
									for (int k2 = 0; k2 < info->DIMENSION; k2++)
									{
										int col = info->Index_Dof[info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
										if (col < 0)
										{
											int col_local = k1 * info->DIMENSION + k2;
											for (int n = 0; n < info->Total_Constraint_to_mesh[Total_mesh]; n++)
											{
												if (info->Controlpoint_of_Element[info->eoi[i][j] * MAX_NO_CP_ON_ELEMENT + k1] == info->Constraint_Node_Dir[n * 2 + 0] && k2 == info->Constraint_Node_Dir[n * 2 + 1])
												{
													double val = -K_linear_coupled[col_local * MAX_KIEL_SIZE + row_local] * info->Value_of_Constraint[n];

													#pragma omp atomic
													info->rhs_vec[row] += val;
												}
											}
										}
									}
							}
						}
				}
			}
		}
	}
}


// 分布荷重の等価節点力を足す
void Add_Equivalent_Nodal_Force_to_F_Vec(information *info)
{
	#pragma omp parallel for
	for (int i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		int index = info->Index_Dof[i];
		if (index >= 0)
			info->rhs_vec[index] += info->Equivalent_Nodal_Force[i];
	}
}